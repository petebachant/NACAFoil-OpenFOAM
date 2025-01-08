#!/usr/bin/env python
"""Script for running the simulation."""

import argparse
import os
import shutil
import subprocess
import sys

import foampy

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "naca_profile",
        default="0012",
        help="Which NACA profile to simulate, e.g., 0012.",
    )
    parser.add_argument(
        "aoa_deg", default=4, help="Angle of attack in degrees."
    )
    parser.add_argument(
        "--overwrite", "-f", action="store_true", default=False
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        default=False,
        help="Do not create a new directory for this case.",
    )
    args = parser.parse_args()
    case_name = f"naca{args.naca_profile}-re2e5-aoa-{args.aoa_deg}"
    if not args.in_place:
        case_dir = os.path.join("cases", case_name)
    else:
        case_dir = "."
    # Copy case files into a case directory, deleting anything that might
    # exist
    # If the case has already been run, we should see it in the results, or
    # maybe we should use DVC to sort this out
    if (
        not args.overwrite
        and not args.in_place
        and os.path.isdir(case_dir)
        and os.listdir(case_dir)
    ):
        print("Case directory is not empty; exiting")
        sys.exit(0)
    if args.overwrite and not args.in_place and os.path.isdir(case_dir):
        # Delete the case and recreate from scratch
        shutil.rmtree(case_dir)
    if not os.path.isdir(case_dir):
        print(f"Creating case directory {case_dir}")
    os.makedirs(case_dir, exist_ok=True)
    system_dir = os.path.join(case_dir, "system")
    os.makedirs(system_dir, exist_ok=True)
    subprocess.run(
        [
            "python",
            "scripts/blockmeshdict.py",
            args.naca_profile,
            args.aoa_deg,
            "--case",
            case_dir,
        ]
    )
    model_names = {
        "k-epsilon": "kEpsilon",
        "laminar": "kEpsilon",
        "k-omega": "kOmega",
    }
    constant_dir = os.path.join(case_dir, "constant")
    os.makedirs(constant_dir, exist_ok=True)
    shutil.copytree("0.org", os.path.join(case_dir, "0"))
    if not args.in_place:
        # All other non template files to copy over
        paths = [
            "constant/transportProperties",
            "constant/turbulenceProperties",
            "system/controlDict",
            "system/createPatchDict",
            "system/fvSchemes",
            "system/fvSolution",
        ]
        for path in paths:
            shutil.copy(path, os.path.join(case_dir, path))
    # Move into the case directory
    print(f"Changing working directory to {case_dir}")
    os.chdir(case_dir)
    # Create the mesh
    foampy.run("blockMesh", overwrite=args.overwrite)
    foampy.run("createPatch", args=["-overwrite"], overwrite=args.overwrite)
    # Run simpleFoam
    foampy.run(
        "simpleFoam",
        overwrite=args.overwrite,
    )
    # Touch case.foam file so we can easily open with ParaView
    subprocess.call(["touch", "case.foam"])
