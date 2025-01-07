#!/usr/bin/env python

from __future__ import division, print_function

import argparse
import os

import numpy as np
from numpy import arctan, cos, linspace, ones, pi, sin, zeros

FRONT_BACK_PATCHES = """
    symFront
    {
        type symmetryPlane;
        faces
        (
            (0 1 5 4)
            (1 2 7 5)
            (7 2 3 8)
            (8 11 10 7)
            (10 9 6 7)
            (6 9 0 4)
        );
    }

    symBack
    {
        type symmetryPlane;
        faces
        (
            (17 13 12 16)
            (19 14 13 17)
            (20 15 14 19)
            (23 20 19 22)
            (22 19 18 21)
            (21 18 16 12)
        );
    }
"""


def gen_blockmeshdict(foil="0012", alpha_deg=4):
    """Write a `blockMeshDict` for a NACA foil at specified angle of attack."""
    # Foil geometry
    c = 1.0  # Geometric chord length
    alpha = np.deg2rad(alpha_deg)  # Angle of attack (in radians)
    NACA = [int(d) for d in foil]  # NACA 4-digit designation
    # Mesh dimensions
    scale = 1  # Scaling factor
    H = 8  # *Half* height of channel
    W = 0.5  # *Half* depth of foil (y-direction)
    D = 16  # Length of downstream section
    # Mesh resolution parameters
    Ni = 400  # Number of interpolation points along the foil
    Nx = 250  # Number of mesh cells along the foil
    ND = 150  # Number of cells in the downstream direction
    NT = 100  # Number of cells the transverse direction
    NW = 1  # Number of cells in the y-direction (along the foil axis)
    # Expansion rates
    ExpT = 500  # Expansion rate in transverse direction
    ExpD = 100  # Expansion rate in the downstream direction
    ExpArc = 50  # Expansion rate along the inlet arc
    # Create a vector with x-coordinates, camber and thickness
    beta = linspace(0, pi, Ni)
    x = c * (0.5 * (1 - cos(beta)))
    z_c = zeros(len(x))
    z_t = zeros(len(x))
    theta = zeros(len(x))
    # Values of m, p and t
    m = NACA[0] / 100
    p = NACA[1] / 10
    t = (NACA[2] * 10 + NACA[3]) / 100
    # Calculate thickness
    # The upper expression will give the airfoil a finite thickness at the trailing
    # edge, witch might cause trouble. The lower expression is corrected to give
    # zero thickness at the trailing edge, but the foil is strictly speaking no
    # longer a proper NACA airfoil.
    #
    # See http://turbmodels.larc.nasa.gov/naca4412sep_val.html
    #     http://en.wikipedia.org/wiki/NACA_airfoil
    # y_t = (t*c/0.2) * (0.2969*(x/c)**0.5 - 0.1260*(x/c) - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1015*(x/c)**4)
    z_t = (t * c / 0.2) * (
        0.2969 * (x / c) ** 0.5
        - 0.1260 * (x / c)
        - 0.3516 * (x / c) ** 2
        + 0.2843 * (x / c) ** 3
        - 0.1036 * (x / c) ** 4
    )
    if p > 0:
        # Calculate camber
        z_c += (m * x / p**2) * (2 * p - x / c) * (x < p * c)
        z_c += (
            (m * (c - x) / (1 - p) ** 2) * (1 + x / c - 2 * p) * (x >= p * c)
        )
        # Calculate theta-value
        theta += arctan((m / p**2) * (2 * p - 2 * x / c)) * (x < p * c)
        theta += arctan((m / (1 - p) ** 2) * (-2 * x / c + 2 * p)) * (
            x >= p * c
        )
    # Calculate coordinates of upper surface
    Xu = x - z_t * sin(theta) - scale / 4.0
    Zu = z_c + z_t * cos(theta)
    # Calculate coordinates of lower surface
    Xl = x + z_t * sin(theta) - scale / 4.0
    Zl = z_c - z_t * cos(theta)
    # Rotate foil to reach specified angle of attack
    upper = np.matrix(
        [[cos(alpha), sin(alpha)], [-sin(alpha), cos(alpha)]]
    ) * np.vstack((Xu, Zu))
    lower = np.matrix(
        [[cos(alpha), sin(alpha)], [-sin(alpha), cos(alpha)]]
    ) * np.vstack((Xl, Zl))
    Xu = upper[0, :].conj().transpose()
    Zu = upper[1, :].conj().transpose()
    Xl = lower[0, :].conj().transpose()
    Zl = lower[1, :].conj().transpose()
    if p > 0:
        # Find index i of max. camber
        C_max_idx = np.where(z_c == max(z_c))[0][0]
    else:
        # Otherwise use location of max. thickness
        C_max_idx = np.where(z_t == max(z_t))[0][0]
    # Move point of mesh "nose"
    NoseX = (-H + Xu[C_max_idx]) * cos(alpha)
    NoseZ = -(-H + Xu[C_max_idx]) * sin(alpha)
    # Calculate the location of the vertices on the positive y-axis and put
    # them in a matrix
    vertices = zeros((12, 3))
    vertices[0, :] = [NoseX[0, 0], W, NoseZ[0, 0]]
    vertices[1, :] = [Xu[C_max_idx, 0], W, H]
    vertices[2, :] = [Xu[-1, 0], W, H]
    vertices[3, :] = [D, W, H]
    vertices[4, :] = [Xu[0, 0], W, Zu[0, 0]]
    vertices[5, :] = [Xu[C_max_idx, 0], W, Zu[C_max_idx, 0]]
    vertices[6, :] = [Xl[C_max_idx, 0], W, Zl[C_max_idx, 0]]
    vertices[7, :] = [Xu[-1, 0], W, Zu[-1, 0]]
    vertices[8, :] = [D, W, Zu[-1, 0]]
    vertices[9, :] = [Xl[C_max_idx, 0], W, -H]
    vertices[10, :] = [Xu[-1, 0], W, -H]
    vertices[11, :] = [D, W, -H]
    # Create vertices for other side (negative y-axis)
    vertices2 = vertices.copy()
    vertices2[:, 1] *= -1
    vertices = np.vstack((vertices, vertices2))
    # Edge 4-5 and 16-17
    pts1 = np.concatenate(
        [
            Xu[1:C_max_idx],
            W * ones(np.shape(Xu[1:C_max_idx])),
            Zu[1:C_max_idx],
        ],
        axis=1,
    )
    pts5 = np.concatenate([pts1[:, 0], -pts1[:, 1], pts1[:, 2]], axis=1)
    # Edge 5-7 and 17-19
    pts2 = np.concatenate(
        [
            Xu[C_max_idx + 1 : Ni - 1],
            W * ones(np.shape(Xu[C_max_idx + 1 : Ni - 1])),
            Zu[C_max_idx + 1 : Ni - 1],
        ],
        axis=1,
    )
    pts6 = np.concatenate([pts2[:, 0], -pts2[:, 1], pts2[:, 2]], axis=1)
    # Edge 4-6 and 16-18
    pts3 = np.concatenate(
        [
            Xl[1:C_max_idx],
            W * ones(np.shape(Xl[1:C_max_idx])),
            Zl[1:C_max_idx],
        ],
        axis=1,
    )
    pts7 = np.concatenate([pts3[:, 0], -pts3[:, 1], pts3[:, 2]], axis=1)
    # Edge 6-7 and 18-19
    pts4 = np.concatenate(
        [
            Xl[C_max_idx + 1 : Ni - 1],
            W * ones(np.shape(Xl[C_max_idx + 1 : Ni - 1])),
            Zl[C_max_idx + 1 : Ni - 1],
        ],
        axis=1,
    )
    pts8 = np.concatenate([pts4[:, 0], -pts4[:, 1], pts4[:, 2]], axis=1)
    # Edge 0-1 and 12-13
    pts9 = np.array([-H * cos(pi / 4) + Xu[C_max_idx, 0], W, H * sin(pi / 4)])
    pts11 = np.array([pts9[0], -pts9[1], pts9[2]])

    # Edge 0-9 and 12-21
    pts10 = np.array(
        [-H * cos(pi / 4) + Xu[C_max_idx, 0], W, -H * sin(pi / 4)]
    )
    pts12 = np.array([pts10[0], -pts10[1], pts10[2]])
    # Calculate number of mesh points along 4-5 and 4-6
    # Nleading = (C_max_idx/Ni)*Nx
    Nleading = int((x[C_max_idx] / c) * Nx)
    # Calculate number of mesh points along 5-7 and 6-7
    Ntrailing = Nx - Nleading
    # Create file content
    txt = "/*--------------------------------*- C++ -*----------------------------------*\\ \n"
    txt += "| =========                 |                                                 | \n"
    txt += "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n"
    txt += "|  \\\\    /   O peration     | Version:  2406                                  | \n"
    txt += "|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      | \n"
    txt += "|    \\\\/     M anipulation  |                                                 | \n"
    txt += "\\*---------------------------------------------------------------------------*/ \n"
    txt += "FoamFile                                                                        \n"
    txt += "{                                                                               \n"
    txt += "    version     2.0;                                                            \n"
    txt += "    format      ascii;                                                          \n"
    txt += "    class       dictionary;                                                     \n"
    txt += "    object      blockMeshDict;                                                  \n"
    txt += "}                                                                               \n"
    txt += "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n"
    txt += "\n"
    txt += "scale %f; \n" % scale
    txt += "\n"
    txt += "vertices \n"
    txt += "( \n"
    for vertex in vertices:
        txt += "    (%f %f %f)\n" % tuple(vertex)
    txt += "); \n"
    txt += "\n"
    txt += "blocks \n"
    txt += "( \n"
    txt += (
        "    hex (4 5 1 0 16 17 13 12)     (%i %i %i) edgeGrading (1 %f %f 1 %f %f %f %f 1 1 1 1) \n"
        % (Nleading, NT, NW, 1 / ExpArc, 1 / ExpArc, ExpT, ExpT, ExpT, ExpT)
    )
    txt += (
        "    hex (5 7 2 1 17 19 14 13)     (%i %i %i) simpleGrading (1 %f 1) \n"
        % (Ntrailing, NT, NW, ExpT)
    )
    txt += (
        "    hex (7 8 3 2 19 20 15 14)     (%i %i %i) simpleGrading (%f %f 1) \n"
        % (ND, NT, NW, ExpD, ExpT)
    )
    txt += (
        "    hex (16 18 21 12 4 6 9 0)     (%i %i %i) edgeGrading (1 %f %f 1 %f %f %f %f 1 1 1 1) \n"
        % (Nleading, NT, NW, 1 / ExpArc, 1 / ExpArc, ExpT, ExpT, ExpT, ExpT)
    )
    txt += (
        "    hex (18 19 22 21 6 7 10 9)    (%i %i %i) simpleGrading (1 %f 1) \n"
        % (Ntrailing, NT, NW, ExpT)
    )
    txt += (
        "    hex (19 20 23 22 7 8 11 10)   (%i %i %i) simpleGrading (%f %f 1) \n"
        % (ND, NT, NW, ExpD, ExpT)
    )
    txt += "); \n"
    txt += "\n"
    txt += "edges \n"
    txt += "( \n"
    txt += "    spline 4 5 \n"
    txt += "        ( \n"
    for pt in np.array(pts1):
        txt += "            (%f %f %f) \n" % tuple(pt)
    txt += "        ) \n"
    txt += "    spline 5 7 \n"
    txt += "        ( \n"
    for pt in np.array(pts2):
        txt += "            (%f %f %f)\n" % tuple(pt)
    txt += "        ) \n"
    txt += "    spline 4 6 \n"
    txt += "        ( \n"
    for pt in np.array(pts3):
        txt += "            (%f %f %f)\n" % tuple(pt)
    txt += "        ) \n"
    txt += "    spline 6 7 \n"
    txt += "        ( \n"
    for pt in np.array(pts4):
        txt += "            (%f %f %f)\n" % tuple(pt)
    txt += "        ) \n"

    txt += "    spline 16 17 \n"
    txt += "        ( \n"
    for pt in np.array(pts5):
        txt += "            (%f %f %f)\n" % tuple(pt)
    txt += "        ) \n"
    txt += "    spline 17 19 \n"
    txt += "        ( \n"
    for pt in np.array(pts6):
        txt += "            (%f %f %f)\n" % tuple(pt)
    txt += "        ) \n"
    txt += "    spline 16 18 \n"
    txt += "        ( \n"
    for pt in np.array(pts7):
        txt += "            (%f %f %f)\n" % tuple(pt)
    txt += "        ) \n"
    txt += "    spline 18 19 \n"
    txt += "        ( \n"
    for pt in np.array(pts8):
        txt += "            (%f %f %f)\n" % tuple(pt)
    txt += "        ) \n"
    txt += "    arc 0 1 (%f %f %f) \n" % tuple(pts9)
    txt += "    arc 0 9 (%f %f %f) \n" % tuple(pts10)
    txt += "    arc 12 13 (%f %f %f) \n" % tuple(pts11)
    txt += "    arc 12 21 (%f %f %f) \n" % tuple(pts12)
    txt += "); \n"
    txt += "\n"
    txt += "boundary \n"
    txt += "( \n"
    txt += "    inlet \n"
    txt += "    { \n"
    txt += "        type patch; \n"
    txt += "        faces \n"
    txt += "        ( \n"
    txt += "            (1 0 12 13) \n"
    txt += "            (0 9 21 12) \n"
    txt += "        ); \n"
    txt += "    } \n"
    txt += "\n"
    txt += "    outlet \n"
    txt += "    { \n"
    txt += "        type patch; \n"
    txt += "        faces \n"
    txt += "        ( \n"
    txt += "            (11 8 20 23) \n"
    txt += "            (8 3 15 20) \n"
    txt += "        ); \n"
    txt += "    } \n"
    txt += FRONT_BACK_PATCHES
    txt += "\n"
    txt += "    topAndBottom \n"
    txt += "    { \n"
    txt += "        type patch; \n"
    txt += "        faces \n"
    txt += "        ( \n"
    txt += "            (3 2 14 15) \n"
    txt += "            (2 1 13 14) \n"
    txt += "            (9 10 22 21) \n"
    txt += "            (10 11 23 22) \n"
    txt += "        ); \n"
    txt += "    } \n"
    txt += "\n"
    txt += "    airfoil \n"
    txt += "    { \n"
    txt += "        type wall; \n"
    txt += "        faces \n"
    txt += "        ( \n"
    txt += "            (5 4 16 17) \n"
    txt += "            (7 5 17 19) \n"
    txt += "            (4 6 18 16) \n"
    txt += "            (6 7 19 18) \n"
    txt += "        ); \n"
    txt += "    } \n"
    txt += "); \n"
    txt += " \n"
    txt += "mergePatchPairs \n"
    txt += "( \n"
    txt += "); \n"
    txt += " \n"
    txt += "// ************************************************************************* // \n"
    # Write to file
    with open("system/blockMeshDict", "w") as f:
        f.write(txt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create blockMeshDict")
    parser.add_argument("foil", help="NACA foil digits")
    parser.add_argument("alpha_deg", help="Angle of attack (deg)")
    args = parser.parse_args()
    print(
        "Generating blockMeshDict for a NACA {} at {} "
        "degrees angle of attack".format(args.foil, args.alpha_deg)
    )
    gen_blockmeshdict(args.foil, float(args.alpha_deg))
