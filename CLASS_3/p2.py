#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
# $Author$
# $Date$
# $URL$
__giturl__ = "$URL$"
import numpy as np

# Example B-spline basis values


N = np.array([0.20, 0.50, 0.30])

# First derivatives
Nx = np.array([-0.8, 0.2, 0.6])
Ny = np.array([-0.5, 0.1, 0.4])

# Second derivatives
Nxx = np.array([1.2, -2.0, 0.8])
Nxy = np.array([0.3, -0.4, 0.1])
Nyy = np.array([0.9, -1.5, 0.6])

# NURBS weights
w = np.array([1.0, 1.0, 1.0])


# Compute Q quantities


Q = w * N
Qx = w * Nx
Qy = w * Ny
Qxx = w * Nxx
Qxy = w * Nxy
Qyy = w * Nyy


# Compute S quantities


S = np.sum(Q)
Sx = np.sum(Qx)
Sy = np.sum(Qy)
Sxx = np.sum(Qxx)
Sxy = np.sum(Qxy)
Syy = np.sum(Qyy)


# Zeroth derivative


R = Q / S


# First derivatives


Rx = (Qx - R * Sx) / S
Ry = (Qy - R * Sy) / S


# Second derivatives


Rxx = (Qxx - 2*Sx*Rx - Sxx*R) / S

Rxy = (Qxy - Sx*Ry - Sy*Rx - Sxy*R) / S

Ryy = (Qyy - 2*Sy*Ry - Syy*R) / S


# output
print("Basis functions (R)")
print(R)

print("\nFirst derivatives")
print("R_x =", Rx)
print("R_y =", Ry)

print("\nSecond derivatives")
print("R_xx =", Rxx)
print("R_xy =", Rxy)
print("R_yy =", Ryy)
