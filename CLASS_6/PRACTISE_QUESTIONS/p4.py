#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
# $Author$
# $Date$
# $URL$
__giturl__ = "$URL$"

import numpy as np
import matplotlib.pyplot as plt
from math import comb

# ==========================================================
# Bernstein Basis
# ==========================================================

def bernstein(p, t):

    B = np.zeros((len(t), p + 1))

    for i in range(p + 1):

        B[:, i] = comb(p, i) * (t ** i) * ((1 - t) ** (p - i))

    return B


# ==========================================================
# Original Cubic Bézier Curve
# ==========================================================

P3 = np.array([
    [0.0, 0.0],
    [1.0, 2.5],
    [3.0, 2.5],
    [4.0, 0.0]
])

t = np.linspace(0, 1, 300)

B3 = bernstein(3, t)

original_curve = B3 @ P3


# ==========================================================
# Method 1
# Simple Order Reduction
# ==========================================================

P2_simple = np.array([
    P3[0],
    (P3[1] + P3[2]) / 2,
    P3[3]
])

B2 = bernstein(2, t)

curve_simple = B2 @ P2_simple


# ==========================================================
# Method 2
# Least Squares Projection
# ==========================================================

P2_LS, residuals, rank, s = np.linalg.lstsq(
    B2,
    original_curve,
    rcond=None
)

curve_LS = B2 @ P2_LS


# ==========================================================
# Compute Errors
# ==========================================================

error_simple = np.sqrt(
    np.mean(np.sum((original_curve - curve_simple) ** 2, axis=1))
)

error_LS = np.sqrt(
    np.mean(np.sum((original_curve - curve_LS) ** 2, axis=1))
)

print("--------------------------------------")
print("Least-Squares Projection Practice")
print("--------------------------------------")

print("\nOriginal Cubic Control Points\n")
print(P3)

print("\nSimple Reduced Quadratic\n")
print(np.round(P2_simple, 4))

print("\nLeast-Squares Quadratic\n")
print(np.round(P2_LS, 4))

print("\nRMS Error")

print("Simple Reduction     =", round(error_simple, 6))
print("Least Squares        =", round(error_LS, 6))


# ==========================================================
# Plot
# ==========================================================

plt.figure(figsize=(9,6))

plt.plot(
    original_curve[:,0],
    original_curve[:,1],
    linewidth=3,
    label="Original Cubic"
)

plt.plot(
    curve_simple[:,0],
    curve_simple[:,1],
    "--",
    linewidth=2,
    label="Simple Order Reduction"
)

plt.plot(
    curve_LS[:,0],
    curve_LS[:,1],
    "-.",
    linewidth=2,
    label="Least-Squares Projection"
)

# Control polygons

plt.plot(
    P3[:,0],
    P3[:,1],
    'o-',
    alpha=0.5
)

plt.plot(
    P2_simple[:,0],
    P2_simple[:,1],
    's--',
    alpha=0.6
)

plt.plot(
    P2_LS[:,0],
    P2_LS[:,1],
    '^-.',
    alpha=0.6
)

plt.title("Least-Squares Projection (IGA Practice)")

plt.axis("equal")

plt.grid(True)

plt.legend()

plt.show()