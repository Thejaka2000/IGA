#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt


#quadratic surface
P = np.array([

    [[0.0, 0.0, 0.0],
     [0.0, 1.0, 0.5],
     [0.0, 2.0, 0.0]],

    [[1.0, 0.0, 1.0],
     [1.0, 1.0, 2.0],
     [1.0, 2.0, 1.0]],

    [[2.0, 0.0, 0.0],
     [2.0, 1.0, 0.5],
     [2.0, 2.0, 0.0]]

])

weights = np.ones((3,3))

# Degrees
p = 2
q = 2

# Knot vectors
Xi  = np.array([0,0,0,1,1,1], dtype=float)
Eta = np.array([0,0,0,1,1,1], dtype=float)

# B-SPLINE BASIS
def bspline_basis(i, p, xi, Xi):

    if p == 0:

        if (Xi[i] <= xi < Xi[i+1]) or \
           (xi == Xi[-1] and Xi[i] <= xi <= Xi[i+1]):
            return 1.0

        return 0.0

    left = 0.0
    right = 0.0

    den1 = Xi[i+p] - Xi[i]
    den2 = Xi[i+p+1] - Xi[i+1]

    if den1 != 0:
        left = ((xi - Xi[i]) / den1) * \
                bspline_basis(i, p-1, xi, Xi)

    if den2 != 0:
        right = ((Xi[i+p+1] - xi) / den2) * \
                 bspline_basis(i+1, p-1, xi, Xi)

    return left + right

# NURBS BASIS
def nurbs_basis(i, j, xi, eta):

    n = P.shape[0]
    m = P.shape[1]

    numerator = (
        bspline_basis(i, p, xi, Xi)
        *
        bspline_basis(j, q, eta, Eta)
        *
        weights[i,j]
    )

    denom = 0.0

    for a in range(n):
        for b in range(m):

            denom += (
                bspline_basis(a,p,xi,Xi)
                *
                bspline_basis(b,q,eta,Eta)
                *
                weights[a,b]
            )

    return numerator / denom

# SURFACE POINT
def surface_point(xi, eta):

    n = P.shape[0]
    m = P.shape[1]

    x = np.zeros(3)

    for i in range(n):
        for j in range(m):

            R = nurbs_basis(i,j,xi,eta)

            x += R * P[i,j]

    return x

# DERIVATIVES
def dx_dxi(xi, eta, h=1e-5):

    return (
        surface_point(xi+h, eta)
        -
        surface_point(xi-h, eta)
    ) / (2*h)

def dx_deta(xi, eta, h=1e-5):

    return (
        surface_point(xi, eta+h)
        -
        surface_point(xi, eta-h)
    ) / (2*h)

def d2x_dxi2(xi, eta, h=1e-5):

    return (
        surface_point(xi+h, eta)
        - 2*surface_point(xi,eta)
        + surface_point(xi-h, eta)
    ) / (h*h)

def d2x_deta2(xi, eta, h=1e-5):

    return (
        surface_point(xi, eta+h)
        - 2*surface_point(xi,eta)
        + surface_point(xi, eta-h)
    ) / (h*h)

def d2x_dxideta(xi, eta, h=1e-5):

    return (
        surface_point(xi+h, eta+h)
        -
        surface_point(xi+h, eta-h)
        -
        surface_point(xi-h, eta+h)
        +
        surface_point(xi-h, eta-h)
    ) / (4*h*h)

# PARAMETRIC LOCATION
xi  = 0.5
eta = 0.5

# FIRST DERIVATIVES
g1 = dx_dxi(xi, eta)
g2 = dx_deta(xi, eta)

# NORMAL VECTOR
cross = np.cross(g1, g2)

n = cross / np.linalg.norm(cross)

# SECOND DERIVATIVES
x11 = d2x_dxi2(xi, eta)
x22 = d2x_deta2(xi, eta)
x12 = d2x_dxideta(xi, eta)

# FIRST FUNDAMENTAL FORM
a11 = np.dot(g1, g1)
a12 = np.dot(g1, g2)
a21 = a12
a22 = np.dot(g2, g2)

A = np.array([
    [a11, a12],
    [a21, a22]
])

# SECOND FUNDAMENTAL FORM
b11 = np.dot(x11, n)
b12 = np.dot(x12, n)
b21 = b12
b22 = np.dot(x22, n)

B = np.array([
    [b11, b12],
    [b21, b22]
])

print("\nFirst Fundamental Form:")
print(A)

print("\nSecond Fundamental Form:")
print(B)

print("\nNormal Vector:")
print(n)

# SURFACE SAMPLING
Nu = 40
Nv = 40

surface_data = []

for u in np.linspace(0,1,Nu):

    for v in np.linspace(0,1,Nv):

        x = surface_point(u,v)

        surface_data.append(x)

    surface_data.append([np.nan,np.nan,np.nan])

surface_data = np.array(surface_data)

np.savetxt(
    "surface.dat",
    surface_data,
    fmt="%.8f"
)

# MATPLOTLIB PLOT
fig = plt.figure(figsize=(8,6))

ax = fig.add_subplot(111, projection='3d')

valid = ~np.isnan(surface_data[:,0])

ax.plot(
    surface_data[valid,0],
    surface_data[valid,1],
    surface_data[valid,2]
)

# Normal vector at parametric point
x0 = surface_point(xi, eta)

ax.quiver(
    x0[0], x0[1], x0[2],
    n[0], n[1], n[2],
    length=0.5
)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.title("NURBS Surface")

plt.show()


# GNUPLOT SCRIPT
gnuplot_script = r'''
set hidden3d
set xlabel "x"
set ylabel "y"
set zlabel "z"

splot "surface.dat" using 1:2:3 with lines lw 1
'''

with open("plot.gnu", "w") as f:
    f.write(gnuplot_script)

print("\nGenerated:")
print("  surface.dat")
print("  plot.gnu")
