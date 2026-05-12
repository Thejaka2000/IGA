#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created: by Thejaka Jayasinghe
import numpy as np
import subprocess
from math import comb

p = 3
dXi = [1, 1, 1, 1]

# Build knot vector
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

# Local coordinate
def local_t(x, e, Xi):
    return (x - Xi[e]) / (Xi[e+1] - Xi[e])

# N^e(t) = 1 - L^e(t) - R^e(t)
def L_e(t):
    return bernstein(3, 0, t)

def R_e(t):
    return bernstein(3, 3, t)

def N_e(t):
    return 1.0 - L_e(t) - R_e(t)

# Output data
xmin = Xi[0]
xmax = Xi[-1]

with open("N_function.dat", "w") as f:
    for x in np.linspace(xmin, xmax, 800):
        N_val = 0.0
        for e in range(len(Xi)-1):
            if Xi[e] <= x <= Xi[e+1] and Xi[e] != Xi[e+1]:
                t = local_t(x, e, Xi)
                N_val += N_e(t)
        f.write(f"{x} {N_val}\n")

with open("R_function.dat", "w") as f:
    for x in np.linspace(xmin, xmax, 800):
        R_val = 0.0
        for e in range(len(Xi)-1):
            if Xi[e] <= x <= Xi[e+1] and Xi[e] != Xi[e+1]:
                t = local_t(x, e, Xi)
                R_val += R_e(t)
        f.write(f"{x} {R_val}\n")

with open("L_function.dat", "w") as f:
    for x in np.linspace(xmin, xmax, 800):
        L_val = 0.0
        for e in range(len(Xi)-1):
            if Xi[e] <= x <= Xi[e+1] and Xi[e] != Xi[e+1]:
                t = local_t(x, e, Xi)
                L_val += L_e(t)
        f.write(f"{x} {L_val}\n")

# GNUplot
caption = "p = 3, knot spans: ΔΞ = {1,1,1,1}"
gp = f"""
set terminal qt size 1000,650
set title "N(Ξ) + R(Ξ) + L(Ξ) using Bezier Extraction"

set xlabel "\\ \\ Ξ"
set ylabel "N(Ξ) + R(Ξ) + L(Ξ)"

set grid
set key outside

set xrange [{xmin}:{xmax}]
set yrange [0:1.05]
set xtics 1
set ytics 0.5

set label 1 "{caption}" at screen 0.5,0.03 center

plot \\
'N_function.dat' using 1:2 with lines lw 3 lc rgb "orange" title "N(Ξ)", \\
'R_function.dat' using 1:2 with lines lw 3 lc rgb "red" title "R(Ξ)", \\
'L_function.dat' using 1:2 with lines lw 3 lc rgb "blue" title "L(Ξ)"

pause -1
"""

subprocess.run(["gnuplot", "-persist"], input=gp, text=True)
