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


# Bernstein basis

def bernstein_basis(p, t):
    B = np.zeros((len(t), p + 1))
    for i in range(p + 1):
        B[:, i] = (
            np.math.comb(p, i)
            * (t ** i)
            * ((1 - t) ** (p - i))
        )
    return B


# Cubic Bézier control points (example)

P3 = np.array([
    [0.0, 0.0],
    [1.0, 2.0],
    [3.0, 2.0],
    [4.0, 0.0]
])

# Sample points on cubic curve
t = np.linspace(0, 1, 200)
B3 = bernstein_basis(3, t)
curve3 = B3 @ P3


# Order reduction
# Cubic -> Quadratic Bézier


B2 = bernstein_basis(2, t)

# Least-squares solution
P2, _, _, _ = np.linalg.lstsq(B2, curve3, rcond=None)

curve2 = B2 @ P2


# Two-element quadratic B-spline approximation
# Split parameter into two elements
t1 = np.linspace(0, 0.5, 100)
t2 = np.linspace(0.5, 1.0, 100)

curve1 = curve3[:100]
curve2_target = curve3[100:]

# Local parameter
u1 = t1 / 0.5
u2 = (t2 - 0.5) / 0.5

B_local = bernstein_basis(2, u1)

# First element
Q1, _, _, _ = np.linalg.lstsq(B_local, curve1, rcond=None)

# Second element
Q2, _, _, _ = np.linalg.lstsq(B_local, curve2_target, rcond=None)

curve_bspline = np.vstack([
    B_local @ Q1,
    B_local @ Q2
])


# Print control points


print("Original cubic control points\n")
print(P3)

print("\nQuadratic Bézier control points\n")
print(np.round(P2,4))

print("\nTwo-element quadratic B-spline control points\n")
print("Element 1")
print(np.round(Q1,4))

print("Element 2")
print(np.round(Q2,4))


# Plot


plt.figure(figsize=(8,6))

plt.plot(curve3[:,0],curve3[:,1],
         label="Original Cubic",
         linewidth=3)

plt.plot(curve2[:,0],curve2[:,1],
         '--',
         label="Quadratic Bézier")

plt.plot(curve_bspline[:,0],curve_bspline[:,1],
         '-.',
         label="Two-element Quadratic B-spline")

plt.scatter(P3[:,0],P3[:,1],s=60)

plt.legend()

plt.axis("equal")
plt.grid(True)

plt.show()
