#!/usr/bin/env python
"""Aggregate all results."""

import os

import numpy as np
import pandas as pd


def load_force_coeffs(casedir: str = ".", timedir: str = "0") -> pd.DataFrame:
    """Load force coefficients time series into a DataFrame."""
    fpath = os.path.join(
        casedir, "postProcessing", "forceCoeffs", timedir, "coefficient.dat"
    )
    data = np.loadtxt(fpath, skiprows=13)
    df = pd.DataFrame()
    df["time"] = data[:, 0]
    df["cl"] = data[:, 3]
    df["cd"] = data[:, 2]
    df["cm"] = data[:, 1]
    return df


def load_yplus(casedir: str = ".", timedir: str = "0") -> pd.DataFrame:
    """Load y+ time series into a DataFrame."""
    fpath = os.path.join(
        casedir, "postProcessing", "yPlus", timedir, "yPlus.dat"
    )
    data = np.loadtxt(fpath, skiprows=2, usecols=[0, 2, 3, 4])
    df = pd.DataFrame()
    df["time"] = data[:, 0]
    df["yplus_min"] = data[:, 1]
    df["yplus_max"] = data[:, 2]
    df["yplus_mean"] = data[:, 3]
    return df


if __name__ == "__main__":
    casedirs = os.listdir("cases")
    rows = []
    for case in casedirs:
        casedir = os.path.join("cases", case)
        profile = case[4:8]
        case_split = case.split("-", maxsplit=3)
        reynolds = float(case_split[1][2:])
        alpha_deg = float(case_split[-1])
        df_force = load_force_coeffs(casedir=casedir)
        df_yplus = load_yplus(casedir=casedir)
        row = dict(
            profile=profile,
            reynolds_number=reynolds,
            alpha_deg=alpha_deg,
            y_plus_mean=df_yplus["yplus_mean"].iloc[-1],
            cl=df_force.iloc[-1]["cl"],
            cd=df_force.iloc[-1]["cd"],
            cm=df_force.iloc[-1]["cm"],
        )
        rows.append(row)
    df = pd.DataFrame(rows)
    print(df)
    df.to_csv("processed/all-simulated.csv", index=False)
