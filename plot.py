#!/usr/bin/env python
"""
This script plots various quantities.
"""

from __future__ import division, print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_foil_perf(foil="0012"):
    df = pd.read_csv("processed/NACA{}.csv".format(foil))
    plt.figure()
    plt.plot(df.alpha_deg, df.cl)

if __name__ == "__main__":
    plot_foil_perf()
    plt.show()