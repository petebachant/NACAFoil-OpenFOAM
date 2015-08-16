#!/usr/bin/env python
"""
This script runs multiple simulations to check sensitivity to parameters.
"""

from __future__ import division, print_function
import numpy as np
import pandas as pd
import os
from subprocess import call

U_infty = 1.0
c = 1.0

def read_force_coeffs():
    """Read force coefficients from output file."""
    data = np.loadtxt("postProcessing/forceCoeffs/0/forceCoeffs.dat",
                      skiprows=9)
    return {"iterations": data.shape[0], "cl": data[-1, 3], 
            "cd": data[-1, 2], "cm": data[-1, 1]}

def read_turbulence_fields():
    """Read sampled turbulence fields."""
    t = max(os.listdir("postProcessing/sets"))
    with open("postProcessing/sets/{}/point_k_omega_epsilon.xy".format(t)) as f:
        line = f.read().split()
    return {"x_turbulence": float(line[0]),
            "k": float(line[1]),
            "omega": float(line[2]),
            "epsilon": float(line[3])}

def set_Re(Re):
    """
    Set Reynolds number (to two significant digits) via kinematic viscosity in 
    the input files.
    """
    print("Setting Reynolds number to {:.1e}".format(Re))
    nu = U_infty*c/Re
    nu = "{:.1e}".format(nu)
    with open("constant/transportProperties.template") as f:
        txt = f.read()
    with open("constant/transportProperties", "w") as f:
        f.write(txt.format(nu=nu))

def read_Re():
    """Read nu from the input files to calculate Reynolds number."""
    with open("constant/transportProperties") as f:
        for line in f.readlines():
            line = line.replace(";", "")
            line = line.split()
            if len(line) > 1 and line[0] == "nu":
                return U_infty*c/float(line[-1])

def alpha_sweep(foil, start, stop, step, Re=2e5):
    """Vary the foil angle of attack and log results."""
    set_Re(Re)
    alpha_list = np.arange(start, stop, step)
    df = pd.DataFrame(columns=["alpha_deg", "cl", "cd", "cm", "Re","iterations",
                               "x_turbulence", "k", "omega", "epsilon"])
    for alpha in alpha_list:
        call("./Allclean")
        call(["./Allrun", foil, str(alpha)])
        d = read_force_coeffs()
        d["alpha_deg"] = alpha
        d.update(read_turbulence_fields())
        d["Re"] = read_Re()
        df = df.append(d, ignore_index=True)
        df.to_csv("processed/NACA{}_{:.1e}.csv".format(foil, Re), index=False)

def Re_alpha_sweep(foil, Re_start, Re_stop, Re_step, alpha_start, alpha_stop,
                   alpha_step):
    """Create a coefficient dataset for a list of Reynolds numbers."""
    pass

if __name__ == "__main__":
    if not os.path.isdir("processed"):
        os.mkdir("processed")
    
    alpha_sweep("0012", 0, 5, 1, Re=2e5)