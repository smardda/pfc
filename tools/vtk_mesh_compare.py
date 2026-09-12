#!/usr/bin/env python3
"""Compare VTK UnstructuredGrid meshes (geometry + topology only).

This tool reads legacy .vtk and .vtu files through the official VTK Python API.
Useful when meshio-based readers fail.

Checks:
    - Number of points
    - Number of cells
    - Cell type distribution
    - Cell connectivity (exact)
    - Point coordinates (tolerant)

Point coordinates are compared using NumPy allclose tolerance test:
    |a - b| <= atol + rtol * |b|
Exit code:
   0 = PASS
   1 = FAIL
"""

# Standard library imports
import argparse

# Third-party imports
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

ATOL = 1e-8
RTOL = 1e-5

IND = "   "


def main() -> int:
    """Run mesh comparison through VTK Python API."""
    args = parse_arguments()
    return compare_meshes(
        args.mesh_a,
        args.mesh_b,
        atol=args.atol,
        rtol=args.rtol,
        skip_points=args.skip_points,
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    descr_msg = (
        "Compare VTK UnstructuredGrid meshes using VTK Python bindings "
        "(more robust than meshio-based readers)"
    )

    ap = argparse.ArgumentParser(
        description=descr_msg,
        epilog=(
            "Point coordinates comparison rule:\n"
            "  |a - b| <= max(rtol * max(|a|, |b|), atol)\n\n"
        ),
    )

    ap.add_argument("mesh_a", help="First mesh file")
    ap.add_argument("mesh_b", help="Second mesh file")

    # Tolerances for coordinates comparison
    ap.add_argument(
        "--atol",
        type=float,
        default=ATOL,
        help=f"absolute tolerance (default = {ATOL})",
    )
    ap.add_argument(
        "--rtol",
        type=float,
        default=RTOL,
        help=f"relative tolerance (default = {RTOL})",
    )

    # Optional: ignore coordinates
    skip_msg = "Skip coordinate comparison (only topology check)."
    ap.add_argument("--skip-points", action="store_true", help=skip_msg)

    return ap.parse_args()


def compare_meshes(
        mesh_a: str,
        mesh_b: str,
        *,
        atol: float = ATOL,
        rtol: float = RTOL,
        skip_points: bool = False,
) -> int:
    """Compare two VTK meshes and return 0 for pass, 1 for fail."""
    vtk_mesh_a = read_mesh(mesh_a)
    vtk_mesh_b = read_mesh(mesh_b)

    # Compare number of mesh points
    if vtk_mesh_a.GetNumberOfPoints() != vtk_mesh_b.GetNumberOfPoints():
        print(
            f"{IND}mesh compare status : ❌ FAILED - #points "
            f"{vtk_mesh_a.GetNumberOfPoints()} vs {vtk_mesh_b.GetNumberOfPoints()}",
        )
        return 1

    # Compare number of mesh cells
    if vtk_mesh_a.GetNumberOfCells() != vtk_mesh_b.GetNumberOfCells():
        print(
            f"{IND}mesh compare status : ❌ FAILED - #cells "
            f"{vtk_mesh_a.GetNumberOfCells()} vs {vtk_mesh_b.GetNumberOfCells()}",
        )
        return 1

    # Convert VTK cell type arrays to NumPy arrays
    ta = vtk_to_numpy(vtk_mesh_a.GetCellTypesArray())
    tb = vtk_to_numpy(vtk_mesh_b.GetCellTypesArray())

    # Compare histograms of cell types
    if not np.array_equal(np.bincount(ta), np.bincount(tb)):
        print(f"{IND}mesh compare status : ❌ FAILED - cell type distribution differs")
        return 1

    # Convert VTK connectivity to Numpy arrays
    ca = vtk_to_numpy(vtk_mesh_a.GetCells().GetData())
    cb = vtk_to_numpy(vtk_mesh_b.GetCells().GetData())

    # Compare connectivity
    if ca.shape != cb.shape or not np.array_equal(ca, cb):
        print(f"{IND}mesh compare status : ❌ FAILED - connectivity differs")
        return 1

    if skip_points:
        # Point coordinates ignored, only topology check
        print(
            f"{IND}mesh compare status : ✅ PASSED - ",
            "mesh topology identical, point coordinates were not compared",
        )
        return 0

    # skip_points is False -> compare point coordinates using tolerances
    pa = vtk_to_numpy(vtk_mesh_a.GetPoints().GetData())
    pb = vtk_to_numpy(vtk_mesh_b.GetPoints().GetData())

    if pa.shape != pb.shape:
        print(f"{IND}mesh compare status : ❌ FAILED - point array shape differs")
        return 1

    # Compute absolute differences and find the max difference
    diff = np.abs(pa - pb)
    maxdiff = float(diff.max()) if diff.size else 0.0

    # Comparison using tolerances:
    # |a - b| <= max(rtol * max(|a|, |b|), atol)
    if not np.allclose(pa, pb, rtol=rtol, atol=atol, equal_nan=True):
        print(
            f"{IND}mesh compare status  : ❌ FAILED - ",
            "point coordinates out of tolerance",
        )
        print(f"{IND}Maximum difference  : {maxdiff:.3e}")
        print(f"{IND}Tolerance            : atol={atol}, rtol={rtol}\n")
        return 1

    print(f"{IND}mesh compare status : ✅ PASSED - point coordinates within tolerance")
    print(f"{IND}Maximum difference  : {maxdiff:.3e}")
    print(f"{IND}Tolerance           : atol={atol}, rtol={rtol}")
    print(f"{IND}Tolerance rule      : |a-b| <= atol + rtol * |b|")

    return 0


def read_mesh(path: str) -> vtk.vtkUnstructuredGrid:
    """Choose the reader depending on file extension."""
    if path.lower().endswith(".vtu"):
        # Reader  for XML format
        r = vtk.vtkXMLUnstructuredGridReader()
    else:
        # Swap reader for Legacy ASCII/binary format
        r = vtk.vtkUnstructuredGridReader()

    r.SetFileName(path)
    r.Update()
    # Return VTK data object (vtkUnstructuredGrid)
    return r.GetOutput()


if __name__ == "__main__":
    raise SystemExit(main())
