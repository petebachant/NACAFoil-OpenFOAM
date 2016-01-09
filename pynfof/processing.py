"""This module contains processing functions."""

import numpy as np
import pandas as pd
import os


def load_force_coeffs_single(timedir):
    """Load force coefficients into a DataFrame for the given time directory."""
    fpath = "postProcessing/forceCoeffs/{}/forceCoeffs.dat".format(timedir)
    data = np.loadtxt(fpath, skiprows=9)
    df = pd.DataFrame()
    df["time"] = data[:, 0]
    df["cl"] = data[:, 3]
    df["cd"] = data[:, 2]
    df["cm"] = data[:, 1]
    return df


def load_force_coeffs(steady=False):
    """Load force coefficients from file.

    If steady, the file from the `0` directory is used, and the last values are
    returned. Otherwise, arrays are loaded from the latest file.
    """
    if steady:
        return load_force_coeffs_single("0")
    else:
        df = pd.DataFrame()
        timedirs = sorted(os.listdir("postProcessing/forceCoeffs"))
        for t in timedirs:
            df = df.append(load_force_coeffs_single(t), ignore_index=True)
    return df.sort_values(by="time")
