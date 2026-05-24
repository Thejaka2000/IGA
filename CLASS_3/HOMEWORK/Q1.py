#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created by Thejaka Jayasinghe

import numpy as np
import matplotlib.pyplot as plt

# Control points
P = np.array([
    [0.0, 0.0],
    [1.0, 2.0],
    [3.0, 3.0],
    [4.0, 0.0]
])

weights = np.array([1.0, 1.0, 1.0, 1.0])

# Quadratic B-spline
p = 2

# Knot vector
Xi = np.array([0, 0, 0, 1, 2, 2, 2], dtype=float)

# Cox-de Boor basis
def bspline_basis(i, p, xi, Xi):

    if p == 0:
        if Xi[i] <= xi < Xi[i+1]:
            return 1.0
        return 0.0
    
    left = 0.0
    right = 0.0

    den1 = Xi[i+p] - Xi[i]
    den2 = Xi[i+p+1] - Xi[i+1]

    if den1 != 0:
        left = ((xi - Xi[i]) / den1) * bspline_basis(i, p-1, xi, Xi)

    if den2 != 0:
        right = ((Xi[i+p+1] - xi) / den2) * bspline_basis(i+1, p-1, xi, Xi)

    return left + right


# NURBS basis
def nurbs_basis(i, p, xi, Xi, w):

    N = np.array([bspline_basis(a, p, xi, Xi)
                  for a in range(len(w))])

    denom = np.sum(N * w)

    if denom == 0:
        return 0.0

    return N[i] * w[i] / denom

# Geometry mapping
def curve_point(xi):

    R = np.array([
        nurbs_basis(i, p, xi, Xi, weights)
        for i in range(len(P))
    ])

    x = np.sum(R[:, None] * P, axis=0)

    return x

# Numerical derivatives
def first_derivative(xi, h=1e-5):

    return (curve_point(xi + h) - curve_point(xi - h)) / (2*h)

def second_derivative(xi, h=1e-5):

    return (
        curve_point(xi + h)
        - 2*curve_point(xi)
        + curve_point(xi - h)
    ) / (h*h)


# Evaluation
xis = np.linspace(0.0, 1.9999, 300)

curve = []
normals = []

for xi in xis:

    x = curve_point(xi)

    g1 = first_derivative(xi)

    norm_g1 = np.linalg.norm(g1)

    t = g1 / norm_g1

    x11 = second_derivative(xi)

    I = np.eye(2)

    curvature_vector = (
        1.0 / norm_g1**2
    ) * ((I - np.outer(t, t)) @ x11)

    kappa = np.linalg.norm(curvature_vector)

    if kappa > 1e-12:
        n = curvature_vector / kappa
    else:
        n = np.array([0.0, 0.0])

    curve.append(x)

    normals.append([
        x[0],
        x[1],
        n[0],
        n[1]
    ])

curve = np.array(curve)
normals = np.array(normals)

# Save curve data
np.savetxt(
    "curve.dat",
    curve,
    header="x y"
)

np.savetxt(
    "normals.dat",
    normals,
    header="x y nx ny"
)

# Plot using matplotlib
plt.figure(figsize=(8,6))

plt.plot(curve[:,0], curve[:,1], linewidth=2)

plt.scatter(P[:,0], P[:,1])

# Normal vectors
skip = 15

for i in range(0, len(normals), skip):

    x, y, nx, ny = normals[i]

    plt.arrow(
        x, y,
        0.2*nx,
        0.2*ny,
        head_width=0.03,
        length_includes_head=True
    )

plt.axis("equal")
plt.xlabel("x")
plt.ylabel("y")

plt.title("NURBS Curve with Normal Vectors")

plt.show()

# GNUplot script
gnu_script = r'''
set size ratio -1
set grid

set xlabel "x"
set ylabel "y"

plot "curve.dat" using 1:2 with lines lw 2 title "NURBS Curve"
'''

with open("plot.gnu", "w") as f:
    f.write(gnu_script)

print("Generated:")
print("  curve.dat")
print("  normals.dat")
print("  plot.gnu")