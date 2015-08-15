#!/usr/bin/env python
"""
This script runs multiple simulations to check sensitivity to parameters.
"""

from __future__ import division, print_function
import numpy as np
import pandas as pd
import os
from subprocess import call


def read_force_coeffs():
    """Read force coefficients from output file."""
    data = np.loadtxt("postProcessing/forceCoeffs/0/forceCoeffs.dat",
                      skiprows=9)
    return {"iterations": data.shape[0], "cl": data[-1, 3], 
            "cd": data[-1, 2], "cm": data[-1, 1]}

def set_nu(nu):
    """Set kinematic viscosity in the input files."""
    pass

def read_Re():
    """Read nu from the input files to calculate Reynolds number."""
    pass

def alpha_sweep(foil, start, stop, step):
    """Vary the foil angle of attack and log results."""
    alpha_list = np.arange(start, stop, step)
    df = pd.DataFrame(columns=["alpha_deg", "cl", "cd", "cm", "iterations"])
    for alpha in alpha_list:
        call("./Allclean")
        call(["./Allrun", foil, str(alpha)])
        d = read_force_coeffs()
        d["alpha_deg"] = alpha
        df = df.append(d, ignore_index=True)
        df.to_csv("processed/NACA{}.csv".format(foil), index=False)

if __name__ == "__main__":
    if not os.path.isdir("processed"):
        os.mkdir("processed")
    
    alpha_sweep("0012", 0, 14, 1)