#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created: by Thejaka Jayasinghe

import numpy as np
import subprocess
from math import comb

p = 3
dXi = [1, 1, 1, 1]

# knot vector
def build_knot_vector(dXi, p):
    knots = [0.0]
    for dx in dXi:
        knots.append(knots[-1] + dx)
    inner = knots[1:-1]
    Xi = [0.0]*(p+1) + inner + [knots[-1]]*(p+1)
    return Xi

Xi = build_knot_vector(dXi, p)

# Bernstein basis
def bernstein(p, k, t):
    return comb(p, k) * (t**k) * ((1 - t)**(p-k))


# Element coordinate
def local_t(x, e, Xi):
    return (x - Xi[e]) / (Xi[e+1] - Xi[e])

# R^e(t) = B_3,3(t) = t^3
def R_e(t):
    return bernstein(3, 3, t)

# Output data
xmin = Xi[0]
xmax = Xi[-1]

with open("R_function.dat", "w") as f:
    for x in np.linspace(xmin, xmax, 800):
        R_val = 0.0
        for e in range(len(Xi)-1):
            if Xi[e] <= x <= Xi[e+1] and Xi[e] != Xi[e+1]:
                t = local_t(x, e, Xi)
                R_val += R_e(t)
        f.write(f"{x} {R_val}\n")

# GNUplot
caption = "Polynomial order: p = 3, knot spans: ΔΞ = {1,1,1,1}"
gp = f"""
set terminal qt size 1000,650
set title "R(Ξ) using Bezier Extraction"

set xlabel "Ξ"
set ylabel "R(Ξ)"

set grid
set key outside

set xrange [{xmin}:{xmax}]
set yrange [0:1.05]
set xtics 1
set ytics 0.5

set label 1 "{caption}" at screen 0.5,0.03 center

plot 'R_function.dat' using 1:2 with lines lw 3 lc rgb "red" title "R(Ξ)"

pause -1
"""

subprocess.run(["gnuplot", "-persist"], input=gp, text=True)