#!/usr/bin/env python3
"""Compare two legacy .vtk files using fieldcompare.

1) Bitwise filecmp check
2) fieldcompare API for field-level comparison with tolerances.

Works ONLY with legacy .vtk (ASCII/BINARY). vtu is NOT supported.
Uses DefaultEquality:
    * exact equality for integer fields
    * fuzzy equality for floating fields using:
        |a - b| <= max(rtol * max(|a|, |b|), atol)

Default Tolerances:
    atol = 1.e-10
    rtol = 1.e-6

Tolerances can be overridden via:
    --atol <value>
    --rtol <value>

Default behaviour:
    - FAIL only if Q or radn fails
    - If Q pass, other fields failures are treated as warnings
Use --strict to require all fields to pass

Exit code: 0 = PASS/WARNING, 1 = FAIL, 2 = usage / file error
"""
# Standard library imports
import argparse
import filecmp
import os

# Third-party imports
import numpy as np
from fieldcompare import FieldComparisonSuite, FieldDataComparator
from fieldcompare.io import read_field_data
from fieldcompare.predicates import DefaultEquality

IND = "   "

ATOL = 1e-10
RTOL = 1e-6

MSG_DIFF_PASS = "filecmp status  : PASSED ✅ Bitwise Identical"
MSG_FC_PASS = "fieldcompare status  : ✅ PASSED - all fields within tolerance"
MSG_FC_FAIL_Q = "fieldcompare status  : ❌ FAILED - Q out of tolerance"
MSG_FC_FAIL_RADN = "fieldcompare status  : ❌ FAILED - radn out of tolerance"
MSG_FC_FAIL_Q_RADN = "fieldcompare status  : ❌ FAILED - Q and radn out of tolerance"
MSG_FC_FAIL_NON_Q = (
    "fieldcompare status  : ❌ FAILED (strict mode) "
    "Non-critical fields are out of tolerance"
)
MSG_FC_WARN_NON_Q = (
    "fieldcompare status  : ⚠️ WARNING - critical fields (Q/radn) within tolerance "
    "but other fields are out"
)


def main() -> int:
    """Run fields comparison through fieldcompare."""
    args = parse_arguments()
    result = compare_files(
        args.file_a,
        args.file_b,
        atol=args.atol,
        rtol=args.rtol,
        strict=args.strict,
    )
    return exit_code(result)


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    descr_msg = (
        "Compare two legacy .vtk files: filecmp first,"
        "then fieldcompare with tolerances."
    )
    parser = argparse.ArgumentParser(
        description=descr_msg,
        epilog=(
            "Comparison rule:\n"
            "  |a - b| <= max(rtol * max(|a|, |b|), atol)\n\n"
            "Requires: fieldcompare (https://gitlab.com/dglaeser/fieldcompare)\n"
            "Recommended installation:\n"
            "   pip install fieldcompare[all]\n\n"
            "Optional (some VTK readers): pip install vtk\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("file_a", help="Reference .vtk file (path can be relative)")
    parser.add_argument("file_b", help="Candidate .vtk file (path can be relative)")
    # Optional tolerances
    parser.add_argument(
        "--atol",
        type=float,
        default=ATOL,
        help=f"absolute tolerance (default = {ATOL})",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=RTOL,
        help=f"relative tolerance (default = {RTOL})",
    )
    # Optional: fail if any field fails, not only Q
    parser.add_argument(
        "--strict",
        action="store_true",
        help="fails if any field fails (default: fails only if Q or radn fails)",
    )

    return parser.parse_args()


def compare_files(
    file_a: str,
    file_b: str,
    *,
    atol: float = ATOL,
    rtol: float = RTOL,
    strict: bool = False,
) -> dict:
    """Compare two VTK files using fieldcompare."""
    print(f"Candidate : {file_b}", flush=True)
    print(f"Reference : {file_a}", flush=True)

    # Check VTK file consistency
    res_check_vfiles = check_vfiles(file_a, file_b)
    if res_check_vfiles is not None:
        return res_check_vfiles

    # First bitwise comparison via filecmp
    bitwise_identical = filecmp.cmp(file_a, file_b, shallow=False)

    if bitwise_identical:
        print(f"{IND}{MSG_DIFF_PASS}", flush=True)
        return make_result("bitwise_passed")

    # In case filecmp fails, use fieldcompare API for field comparison with tolerances
    tolerance_msg = (
        f"{IND}Tolerance            : atol={atol}, rtol={rtol}\n"
        f"{IND}Rule                 : |a-b| <= max(rtol * max(|a|, |b|), atol)"
    )

    try:
        predicate = DefaultEquality(abs_tol=atol, rel_tol=rtol)

        max_diffs = {}
        tol_bins = {}

        def selector(source_field, reference_field):
            a = np.asarray(source_field.values)
            b = np.asarray(reference_field.values)

            # Same tolerance rule used by fieldcompare
            tolerance = np.maximum(rtol*np.maximum(np.abs(a), np.abs(b)), atol)

            # Absolute difference
            diff = np.abs(a-b)

            # Maximum absolute difference
            idx = np.unravel_index(np.argmax(diff),diff.shape)

            # Compare differences with the allowed tolerance
            ratio = diff/tolerance
            idx_ratio = np.unravel_index(np.argmax(ratio),ratio.shape)

            # Values which fail the tolerance test
            failed = ratio > 1.0

            max_diffs[source_field.name] = {
                "max_abs_diff_cell": int(idx[0]),
                "candidate_value": source_field.values[idx],
                "reference_value": reference_field.values[idx],
                "max_abs_diff": diff[idx],
                "max_tol_violation_cell": int(idx_ratio[0]),
                "max_tol_ratio_candidate": source_field.values[idx_ratio],
                "max_tol_ratio_reference": reference_field.values[idx_ratio],
                "max_tol_ratio": ratio[idx_ratio],
                "failed_count": np.count_nonzero(failed),
                "total_cell_count": diff.shape[0],
            }

            tol_bins[source_field.name] = {
                "1_2": np.count_nonzero((ratio>1.0) & (ratio <= 2.0) ),
                "2_10": np.count_nonzero((ratio>2.0) & (ratio <= 10.0) ),
                "10_100": np.count_nonzero((ratio>10.0) & (ratio<=100)),
                "gt_100": np.count_nonzero(ratio>100),
            }

            return predicate

        fields_a = read_field_data(file_a)
        fields_b = read_field_data(file_b)

        # FieldDataComparator args in order are: source, reference
        comparator = FieldDataComparator(fields_b, fields_a)

        result = comparator(
            predicate_selector=selector,
            fieldcomp_callback=lambda _: None,
        )

    except (SystemExit, Exception):  # noqa: BLE001
        print(f"{IND}fieldcompare: ❌ ERROR - could not read or compare file",
              flush=True,
        )
        return make_result("fieldcompare_error")

    # The result is a suite of field comparison that we can loop over
    if result:  # Status PASS
        print(f"{IND}{MSG_FC_PASS}", flush=True)
        print(f"{IND}fieldcompare report  : {result.report}", flush=True)
        print(f"{tolerance_msg}",flush=True)

        return make_result(
            status="fieldcompare_passed",
            failed_fields=len(result.failed),
            passed_fields=len(result.passed),
            report=result.report,
        )

    # There are failures
    if "Domain equality check failed" in result.report:
        print(f"{IND}fieldcompare: ❌ ERROR - could not read or compare files",
              flush=True)
        return make_result(status="fieldcompare_error")

    q_failed = any(comp.name.startswith("Q @") for comp in result.failed)
    radn_failed = any(comp.name.startswith("radn @") for comp in result.failed)

    if q_failed and radn_failed:
        # Q field  out of tolerance
        err_msg =  MSG_FC_FAIL_Q_RADN
        status = "q_radn_failed"

    elif q_failed:
        # Q field  out of tolerance
        err_msg = MSG_FC_FAIL_Q
        status = "q_failed"

    elif radn_failed:
        # radiation field out of tolerance
        err_msg = MSG_FC_FAIL_RADN
        status = "radn_failed"

    # If we arrive here, Q passed and only non-Q fields failed.
    elif strict:
        # Q and radn passed but others fields out of tolerances
        err_msg = MSG_FC_FAIL_NON_Q
        status = "strict_failed"

    else:
        # Q and radn passed but others fields out of tolerances
        # In non-strict mode this is only a warning
        err_msg = MSG_FC_WARN_NON_Q
        status = "warning"

    print_fields_failed(err_msg, tolerance_msg, max_diffs, tol_bins, result)

    return make_result(
        status=status,
        failed_fields=len(result.failed),
        passed_fields=len(result.passed),
        report=result.report,
    )


# helper functions
def print_fields_failed(
        err_msg: str,
        tolerance_msg: str,
        max_diffs: dict,
        tol_bins: dict,
        result: FieldComparisonSuite
) -> None:
    """Print failed fields."""

    print()
    print(f"{IND}{err_msg}", flush=True)
    print(f"{IND}fieldcompare report  : {result.report}", flush=True)
    print(f"{tolerance_msg}", flush=True)
    print()
    print(f"{IND}Diagnostics for failed fields", flush=True)

    print(
        f"{IND}{'Field':<10} "
        f"{'Failed':>9}"
        f"{'MaxΔ-Cell':>12} "
        f"{'max|Δ|':>12} "
        f"{'Candidate':>14} "
        f"{'Reference':>14} "
        f"{'Max-Δ/Tol-Cell':>15} "
        f"{'xTol':>11} "
        f"{'Candidate':>14} "
        f"{'Reference':>14} ",
        flush=True
    )

    print(f"{IND}{'_' * 133}", flush=True)

    for field_name, data in max_diffs.items():
        if not any(comp.name.startswith(f"{field_name} @") for comp in result.failed):
            continue

        failed_str = f"{data['failed_count']}/{data['total_cell_count']}"
        print(
            f"{IND}{field_name:<10} "
            f"{failed_str:>9}  "
            f"{data['max_abs_diff_cell']:>10d} "
            f"{data['max_abs_diff']:>12.6g} "
            f"{data['candidate_value']:>14.6g} "
            f"{data['reference_value']:>14.6g} "
            f"{data['max_tol_violation_cell']:>15d} "
            f"{data['max_tol_ratio']:>11.3g} "
            f"{data['max_tol_ratio_candidate']:>14.6g} "
            f"{data['max_tol_ratio_reference']:>14.6g} ",
            flush=True,
        )

    print()
    print(
        f"{IND}Failed-value distribution by tolerance ratio "
        f"(xTol = |Δ| / tolerance)"
    )
    for field_name, bins in tol_bins.items():
        if not any(comp.name.startswith(f"{field_name} @") for comp in result.failed):
            continue
        #if field_name in ("Q", "radn") :
        print(
            f"{IND}{field_name:<8}    "
            f"1<xTol<=2:  {bins['1_2']}     "
            f"2<xTol<=10:  {bins['2_10']}     "
            f"10<xTol<=100:  {bins['10_100']}     "
            f"xTol>100:  {bins['gt_100']} ",
            flush=True,
        )

    print()
    print(f"{IND}fieldcompare first out-of-tolerance values",flush=True)
    for comp in result.failed:
        print(f"{IND * 2}[FAILED]: {comp.name}", flush=True)
        print(f"{IND * 5} {comp.report}", flush=True)


def check_vfiles(file_a: str, file_b: str) -> dict | None:
    """Return an error result if files are invalid, otherwise None."""
    if not os.path.exists(file_a):
        print("❌ ERROR reference file not found:", file_a, flush=True)
        return make_result(status="file_a_not_found")

    if not os.path.exists(file_b):
        print("❌ ERROR candidate file not found:", file_b, flush=True)
        return make_result(status="file_b_not_found")

    if not file_a.lower().endswith(".vtk"):
        print(" ❌ ERROR file ", file_a, "is not a legacy .vtk", flush=True)
        return make_result(status="file_a_not_legacy")

    if not file_b.lower().endswith(".vtk"):
        print("❌ ERROR file ", file_b, "is not a legacy .vtk", flush=True)
        return make_result(status="file_b_not_legacy")

    return None


def make_result(
    status: str = "",
    failed_fields: int = 0,
    passed_fields: int = 0,
    report: str = "",
) -> dict:
    """Compose the returning result from vtk_fieldcompare."""
    return {
        "status": status,
        "failed_fields": failed_fields,
        "passed_fields": passed_fields,
        "report": report,
    }


def exit_code(result: dict) -> int:
    """Compute the int exit code from result dictionary."""
    if result["status"] in {
        "file_a_not_found",
        "file_b_not_found",
        "file_a_not_legacy",
        "file_b_not_legacy",
        "fieldcompare_error",
    }:
        return 2
    if result["status"] in {
        "q_radn_failed",
        "q_failed",
        "radn_failed",
        "strict_failed",
    }:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


