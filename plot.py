#!/usr/bin/env python
"""
This script plots various quantities.
"""

from __future__ import division, print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

ylabels = {"cl": r"$C_l$", "cd": r"$C_d$", "cl/cd": r"$C_l/C_d$", "k": "$k$",
           "omega": r"$\omega$", "epsilon": r"$\epsilon$"}

def plot_foil_perf(quantity="cl/cd", foil="0012", Re=2e5):
    df = pd.read_csv("processed/NACA{}_{:.1e}.csv".format(foil, Re))
    plt.figure()
    if quantity == "cl/cd":
        q = df.cl/df.cd
    else:
        q = df[quantity]
    plt.plot(df.alpha_deg, q, "-o")
    plt.xlabel(r"$\alpha$ (deg)")
    plt.ylabel(ylabels[quantity])
    plt.grid(True)
    plt.tight_layout()

if __name__ == "__main__":
    try:
        import seaborn
        seaborn.set(style="white", context="notebook", font_scale=1.5)
    except ImportError:
        print("Could not import seaborn for plot styling. Try")
        print("\n    conda install seaborn\n\nor")
        print("\n    pip install seaborn\n")
        
    parser = argparse.ArgumentParser(description="Plotting results")
    parser.add_argument("quantity", nargs="?", default="cl/cd", 
                        help="Which quantity to plot", 
                        choices=["cl", "cd", "cl/cd", "k", "omega", "epsilon"])
    parser.add_argument("--foil", "-f", help="Foil", default="0012")
    parser.add_argument("--Reynolds", "-R", help="Reynolds number", default=2e5)
    parser.add_argument("--save", "-s", action="store_true", help="Save plots")
    parser.add_argument("--noshow", action="store_true", default=False, 
                        help="Do not show")
    args = parser.parse_args()

    plot_foil_perf(args.quantity, args.foil, args.Reynolds)
    
    if args.save:
        if not os.path.isdir("figures"):
            os.mkdir("figures")
        plt.savefig("figures/{}.pdf".format(args.quantity))
        
    if not args.noshow:
        plt.show()
