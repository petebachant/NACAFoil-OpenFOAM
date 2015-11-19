#!/usr/bin/env python
"""
This script plots various quantities.
"""

from __future__ import division, print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import pynfof.processing as proc

labels = {"cl": r"$C_l$", "cd": r"$C_d$", "cl/cd": r"$C_l/C_d$", "k": "$k$",
          "omega": r"$\omega$", "epsilon": r"$\epsilon$",
          "alpha_deg": r"$\alpha$ (deg)", "cm": r"$C_m$", "cn": "$C_n$",
          "cc": "$C_c$"}


def plot_time_series(quantity="cl"):
    """
    Plot specified quantity over time. Can be used to visualize convergence.
    """
    df = proc.load_force_coeffs()
    if quantity == "cl/cd":
        q = df.cl/df.cd
    else:
        q = df[quantity]
    plt.figure()
    plt.plot(df.time[5:], q[5:])
    plt.xlabel(r"$t$")
    plt.ylabel(labels[quantity])
    plt.grid(True)
    plt.tight_layout()


def plot_foil_perf(quantity="cl/cd", foil="0012", Re=2e5, x="alpha_deg",
                   ax=None, marker="-o"):
    df = pd.read_csv("processed/NACA{}_{:.1e}.csv".format(foil, Re))
    alpha = np.deg2rad(df.alpha_deg)
    df["cn"] = df.cl*np.cos(alpha) - df.cd*np.sin(alpha)
    df["cc"] = df.cl*np.sin(alpha) - df.cd*np.cos(alpha)
    if ax is None:
        fig, ax = plt.subplots()
    if quantity == "cl/cd":
        q = df.cl/df.cd
    else:
        q = df[quantity]
    ax.plot(df[x], q, "-o", label="NACA " + foil)
    ax.set_xlabel(labels[x])
    ax.set_ylabel(labels[quantity])
    ax.grid(True)
    try:
        fig.tight_layout()
    except UnboundLocalError:
        pass


def plot_multiple_foils(quantity="cl/cd", foils=["0012", "0021"], Re=2e5,
                        x="alpha_deg", save=False):
    """Plot performance for multiple foils."""
    fig, ax = subplots()
    for foil in foils:
        plot_foil_perf(quantity=quantity, foil=foil, Re=Re, x=x, ax=ax)
    ax.legend(loc="best")
    fig.tight_layout()
    if save:
        fig.savefig("figures/NACA_" + "-".join(foils) + ".pdf")
        fig.savefig("figures/NACA_" + "-".join(foils) + ".png", dpi=300)


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
                        choices=["cl", "cd", "cm", "cl/cd", "k", "omega",
                                 "epsilon", "cn", "cc"])
    parser.add_argument("-x", help="Quantity on x-axis", default="alpha_deg")
    parser.add_argument("--foil", "-f", help="Foil", default="0012")
    parser.add_argument("--Reynolds", "-R", help="Reynolds number", default=2e5)
    parser.add_argument("--save", "-s", action="store_true", help="Save plots")
    parser.add_argument("--noshow", action="store_true", default=False,
                        help="Do not show")
    parser.add_argument("--timeseries", "-t", action="store_true",
                        default=False, help="Plot time series data")
    args = parser.parse_args()

    if args.timeseries:
        plot_time_series(args.quantity)
    else:
        plot_foil_perf(args.quantity, args.foil, float(args.Reynolds),
                       x=args.x)

    if args.save:
        if not os.path.isdir("figures"):
            os.mkdir("figures")
        plt.savefig("figures/{}.pdf".format(args.quantity))

    if not args.noshow:
        plt.show()
