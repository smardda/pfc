#!/usr/bin/env python3
# Developed by: Francesca Anna Riccio
"""vtk_folder_compare.py - Compare VTK files between two folders.

By default, this script compares every .vtk file found under DIR_A with the file at the
same relative path under DIR_B.
Use --files to compare a selected substet of VTK files from DIR_A.

Comparison strategy:

1) Bitwise diff (-q) for exact matches
2) For files under an 'H/' folder:
      use compare_meshes() from vtk_mesh_compare.py
3) For all other differing VTK files:
      use fieldcompare API for field-level comparison with tolerances.

Default regression rule:
    - FAIL only if Q or radn fails
    - Other field failures are treated as warnings
Use --strict to require all fields to pass.

Exit code: 0 = PASS/WARNING, 1 = FAIL, 2 = usage/file error
"""

# Standard library imports
import argparse
import filecmp
import os

# Local imports
import vtk_fieldcompare
import vtk_mesh_compare

FIELD_ATOL = 1e-8
FIELD_RTOL = 1e-5

MESH_ATOL = 1e-8
MESH_RTOL = 1e-5

IND = "   "
MSG_FILECMP_NOPASS = "filecmp ⚠️ files not bitwise identical"


def main() -> int:
     """Run folder comparison of VTK files."""
     args = parse_arguments()
     return compare_folders(
          args.dir_a,
          args.dir_b,
          files=args.files,
          atol=args.atol,
          rtol=args.rtol,
          strict=args.strict,
     )


def parse_arguments()  -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare .vtk files recursively between two folders: "
            "diff -q first, then mesh compare or fieldcompare with tolerances."
        ),
        epilog=(
            "Field comparison rule:\n"
            "   |a - b| <= max(rtol * max(|a|, |b|), atol)\n\n"
            "Requirements:\n"
            "   - fieldcompare (https://gitlab.com/dglaeser/fieldcompare)\n"
            "       pip install fieldcompare[all]\n"
            "   - vtk_mesh_compare.py must be importable from this script location\n"
            "   - VTK Python bindings: pip install vtk\n\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("dir_a", help="Reference folder")
    parser.add_argument("dir_b", help="Candidate results folder")

    parser.add_argument(
        "--files",
        nargs="+",
        help="subset of VTK files in dir_a to compare",
    )

    parser.add_argument(
        "--atol",
        type=float,
        default=FIELD_ATOL,
        help=f"field absolute tolerance (default = {FIELD_ATOL})",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=FIELD_RTOL,
        help=f"field relative tolerance (default = {FIELD_RTOL})",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="fails if any field fails (default: fails only if Q or radn fails)",
    )

    return  parser.parse_args()

def compare_folders(
    dir_a: str,
    dir_b: str,
    *,
    files: list[str] | None = None,
    atol: float = FIELD_ATOL,
    rtol: float = FIELD_RTOL,
    strict: bool = False,
) -> int:
    """Compare selected VTK files in two regression folders."""
    dir_a = os.path.abspath(dir_a)
    dir_b = os.path.abspath(dir_b)
    print(f"Candidate results folder : {dir_b}")
    print(f"Reference results folder : {dir_a}")

    if not os.path.isdir(dir_a):
        print(f"ERROR: DIR_A doesn't exist: {dir_a}")
        return 2

    if not os.path.isdir(dir_b):
        print(f"ERROR: DIR_B doesn't exist: {dir_b}")
        return 2

    atol_mesh = max(atol, MESH_ATOL)
    rtol_mesh = max(rtol, MESH_RTOL)

    # Counters
    compared = 0
    identical = 0
    acceptable_mesh = 0
    acceptable_field = 0
    warning_field = 0
    fieldcompare_errors = 0
    failed = 0
    failed_field = 0
    failed_mesh = 0
    missing_files = 0

    invalid_files = []
    if files is not None:
        for fname in files:
            fullpath_a = os.path.join(dir_a, fname)
            if not os.path.isfile(fullpath_a):
                invalid_files.append(fname)

    if invalid_files:
        for fname in invalid_files:
            print(
                f"ERROR: file {fname} requested with --files "
                f"is not a file in DIR_A: {dir_a}",
            )
        return 2

    if files is not None:
        vtk_files = files
    else:
        vtk_files = []
        # walk through the VTK files in subfolders in dir_a
        # for root, dirs, files in os.walk(dir_a):

        for root, _, fnames in os.walk(dir_a):
            for fname in fnames:
                if not fname.endswith(".vtk"):
                    continue

                fullpath_a = os.path.join(root, fname)
                relpath = os.path.relpath(fullpath_a, dir_a)

                vtk_files.append(relpath)

    for relpath in vtk_files:
        fullpath_a = os.path.join(dir_a, relpath)
        fullpath_b = os.path.join(dir_b, relpath)
        if os.path.islink(fullpath_a):
            print(relpath, " Skipping symbolic link ", fullpath_a, flush=True)
            continue
        if os.path.islink(fullpath_b):
            print(relpath, " Skipping symbolic link ", fullpath_b, flush=True)
            continue

        if not os.path.exists(fullpath_b):
            missing_files += 1
            print(relpath, " Missing file in dir_b", fullpath_b, flush=True)
            continue

        #  Bitwise comparison (fast check)
        bitwise_identical = filecmp.cmp(fullpath_a, fullpath_b, shallow=False)

        compared += 1
        if bitwise_identical:
            identical += 1
            continue

        print()
        print(relpath, MSG_FILECMP_NOPASS, flush=True)

        # Fallback for file differ
        # H/*.vtk files are mesh/domain files, use vtk_mesh_compare with
        # tolerances, based on VTK Python API, meshio can fail in this case.
        # #If files differ and are under an H folder, use mesh comparison
        if os.path.basename(os.path.dirname(fullpath_a)) == "H":
            print(" -> using vtk_mesh_compare", flush=True)
            mesh_status = vtk_mesh_compare.compare_meshes(
                fullpath_a,
                fullpath_b,
                atol=atol_mesh,
                rtol=rtol_mesh,
                skip_points = False,
            )

            if mesh_status == 0:
                acceptable_mesh += 1
            else:
                failed_mesh += 1
            continue

        # Fallback for file differ
        # S/*.vtk G/*.vtk P/*.vtk files contain fields.
        # Use fieldcompare for field-level comparison with tolerances
        print(" -> using vtk_fieldcompare", flush=True)
        r_fieldcompare = vtk_fieldcompare.compare_files(
            fullpath_a,
            fullpath_b,
            atol=atol,
            rtol=rtol,
            strict=strict,
        )

        status = r_fieldcompare["status"]
        if status == "bitwise_passed":
            identical += 1
        elif status == "fieldcompare_passed":
            acceptable_field += 1
        elif status in {"q_failed", "radn_failed", "q_radn_failed", "strict_failed"}:
            failed_field += 1
        elif status == "warning":
            warning_field += 1
        elif status == "fieldcompare_error":
            fieldcompare_errors += 1
            print("   Trying mesh_compare fallback...", flush=True)
            mesh_status = vtk_mesh_compare.compare_meshes(
                fullpath_a,
                fullpath_b,
                atol=atol_mesh,
                rtol=rtol_mesh,
                skip_points=False,
            )
            if mesh_status == 0:
                acceptable_mesh += 1
            else:
                failed_mesh += 1
        else:
            failed_field += 1

    selected_files = compared + missing_files
    count_figures = len(str(selected_files))
    print("\n============== VTK COMPARISON REPORT ==================")
    print(f"Candidate results folder   : {dir_b}")
    print(f"Reference results folder   : {dir_a}")

    print(f"{IND}Selected .vtk files in reference folder    : {selected_files:>{count_figures}}")
    print(f"{IND}Total number of .vtk files compared        : {compared:>{count_figures}}")
    print(f"{IND}✅ Bitwise identical                       : {identical:>{count_figures}}")

    if acceptable_mesh:
        print(f"{IND}✅ Passed via mesh comparison              : {acceptable_mesh:>{count_figures}}")

    if acceptable_field:
        print(
            f"{IND}✅ Passed via field comparison         : {acceptable_field:>{count_figures}}"
            "  ( all fields passed within tolerance)",
        )

    if warning_field:
        print(
            f"{IND}⚠️ Passed with warnings                   : {warning_field:>{count_figures}}"
            "  (Q/radn are within tolerance, other fields are out)",
        )

    failed += failed_mesh + failed_field
    print(f"{IND}❌ Comparison failures                     : {failed:>{count_figures}}")

    if failed_field:
        print(f"{IND * 2}of which failed via field comparison    : {failed_field:>{count_figures}}")

    if failed_mesh:
        print(f"{IND * 2}of which failed via mesh comparison    : {failed_mesh:>{count_figures}}")

    if missing_files:
        print(f"{IND }❌ Missing in candidate results folder     : {missing_files:>{count_figures}}")

    if fieldcompare_errors:
        print(
            f"{IND * 2}fieldcompare errors with mesh fallback : "
            f"{fieldcompare_errors:>{count_figures}}",
        )

    if compared > 0 and identical == compared:
        print(f"{IND}All compared files were bitwise identical")
        print(f"{IND}Field tolerances were not used")

    if failed_field or acceptable_field or warning_field:
        print(f"{IND}Field tolerances were used:")
        print(f"{IND * 2}atol = {atol}")
        print(f"{IND * 2}rtol = {rtol}")
        print(f"{IND * 2}rule : |a-b| <= max(rtol * max(|a|, |b|), atol)")

    if failed_mesh or acceptable_mesh:
        print(f"{IND}Mesh tolerances were used:")
        print(f"{IND * 2}atol = {atol_mesh}")
        print(f"{IND * 2}rtol = {rtol_mesh}")
        print(f"{IND * 2}rule : |a-b| <= atol + rtol * |a|")

    return 0 if failed == 0 and missing_files == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
