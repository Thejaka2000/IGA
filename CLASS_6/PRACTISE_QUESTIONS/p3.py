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


# Bernstein basis


def bernstein_basis(p, t):
    t = np.asarray(t)
    B = np.zeros((len(t), p + 1))

    for i in range(p + 1):
        B[:, i] = comb(p, i) * (t**i) * (1 - t)**(p - i)

    return B



# Original quadratic B-spline
# (represented by two quadratic Bézier elements)


Q1 = np.array([
    [0.0, 0.0],
    [1.0, 2.0],
    [2.0, 1.5]
])

Q2 = np.array([
    [2.0, 1.5],
    [3.0, 1.0],
    [4.0, 0.0]
])

u = np.linspace(0, 1, 100)

B2 = bernstein_basis(2, u)

curve1 = B2 @ Q1
curve2 = B2 @ Q2

original = np.vstack((curve1, curve2))


# Knot removal
# Approximate by ONE quadratic Bézier


u_all = np.linspace(0, 1, len(original))

Bfit = bernstein_basis(2, u_all)

P_quad, _, _, _ = np.linalg.lstsq(Bfit, original, rcond=None)

quad_curve = Bfit @ P_quad


# Order Elevation
# Quadratic -> Cubic


P_cubic = np.zeros((4, 2))

P_cubic[0] = P_quad[0]

P_cubic[1] = (1/3) * P_quad[0] + (2/3) * P_quad[1]

P_cubic[2] = (2/3) * P_quad[1] + (1/3) * P_quad[2]

P_cubic[3] = P_quad[2]


# Remove another knot
# (one cubic Bézier element)


B3 = bernstein_basis(3, u_all)

cubic_curve = B3 @ P_cubic


# Print results


print("\nOriginal quadratic B-spline\n")
print("Element 1")
print(Q1)

print("\nElement 2")
print(Q2)

print("\nQuadratic after knot removal\n")
print(np.round(P_quad,4))

print("\nCubic after order elevation\n")
print(np.round(P_cubic,4))


# Plot


plt.figure(figsize=(9,6))

plt.plot(original[:,0],
         original[:,1],
         linewidth=3,
         label="Original quadratic B-spline (2 elements)")

plt.plot(quad_curve[:,0],
         quad_curve[:,1],
         "--",
         linewidth=2,
         label="After removing one knot")

plt.plot(cubic_curve[:,0],
         cubic_curve[:,1],
         "-.",
         linewidth=2,
         label="Order elevated (p=3) + remove two knots")

plt.scatter(Q1[:,0],Q1[:,1],s=60)
plt.scatter(Q2[:,0],Q2[:,1],s=60)

plt.scatter(P_quad[:,0],P_quad[:,1],
            marker='s',
            s=60)

plt.scatter(P_cubic[:,0],P_cubic[:,1],
            marker='^',
            s=70)

plt.axis("equal")
plt.grid(True)
plt.legend()
plt.title("Knot Removal")

plt.show()
