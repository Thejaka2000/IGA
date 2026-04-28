#!/usr/bin/env python3
import numpy as np
import subprocess

# USER INPUT
# Change knot spans here
# Example:
# ΔΞ = {1,1,1}
# ΔΞ = {0.5,1.5,1}
# ΔΞ = {2,1,3}

dXi = [1, 1, 1]

p = 2   # quadratic


# Build open knot vector from knot spans
def build_knot_vector(dXi, p):

    knots = [0.0]

    for dx in dXi:
        knots.append(knots[-1] + dx)

    # internal breakpoints
    inner = knots[1:-1]

    Xi = [0.0]*(p+1) + inner + [knots[-1]]*(p+1)

    return Xi


Xi = build_knot_vector(dXi, p)

# Example for [1,1,1]
# Xi = [0,0,0,1,2,3,3,3]


# Cox-de Boor recursion
def N(i, p, x, Xi):

    if p == 0:
        if (Xi[i] <= x < Xi[i+1]) or (x == Xi[-1] and Xi[i] <= x <= Xi[i+1]):
            return 1.0
        return 0.0

    left = 0.0
    right = 0.0

    d1 = Xi[i+p] - Xi[i]
    d2 = Xi[i+p+1] - Xi[i+1]

    if d1 > 0:
        left = (x - Xi[i]) / d1 * N(i, p-1, x, Xi)

    if d2 > 0:
        right = (Xi[i+p+1] - x) / d2 * N(i+1, p-1, x, Xi)

    return left + right


# Number of basis functions
n_basis = len(Xi) - p - 1


# Write data file
xmin = Xi[0]
xmax = Xi[-1]

with open("bspline.dat", "w") as f:

    for x in np.linspace(xmin, xmax, 600):

        vals = [x]

        for i in range(n_basis):
            vals.append(N(i, p, x, Xi))

        f.write(" ".join(map(str, vals)) + "\n")


# GNUplot
span_text = ", ".join(str(v) for v in dXi)
caption = f"Polynomial order: p = {p}, knot spans: ΔΞ = {{{span_text}}}"
gp = f"""
set terminal qt size 900,600
set title "Quadratic B-Spline Basis Functions"
set xlabel "Ξ"
set ylabel "N"
set grid
set key outside
set xrange [{xmin}:{xmax}]
set yrange [0:1.05]

set label 1 "{caption}" at screen 0.5,0.03 center

plot for [i=2:{n_basis+1}] 'bspline.dat' using 1:i with lines lw 2 title sprintf("N_%d,2", i-2)

pause -1
"""

subprocess.run(["gnuplot", "-persist"], input=gp, text=True)
