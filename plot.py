#!/usr/bin/env python
"""
This script plots various quantities.
"""

from __future__ import division, print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_foil_perf(foil="0012", Re=2e5):
    df = pd.read_csv("processed/NACA{}_{:.1e}.csv".format(foil, Re))
    plt.figure()
    plt.plot(df.alpha_deg, df.cl/df.cd, "-o")
    plt.xlabel(r"$\alpha$ (deg)")
    plt.ylabel(r"$C_l/C_d$")
    plt.tight_layout()

if __name__ == "__main__":
    plot_foil_perf()
    plt.show()
