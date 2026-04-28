#!/usr/bin/env python3
import numpy as np
import subprocess

# =====================================================
# Cubic B-spline basis functions
# p = 3
# Knot spans: ΔΞ = {1,2,3,1}
# Show x-axis only at knot locations: 1,3,6,7
# Vertical dotted lines at knots
# =====================================================

p = 3
dXi = [1, 2, 3, 1]


# =====================================================
# Build open knot vector
# =====================================================
def build_knot_vector(dXi, p):

    knots = [0.0]

    for dx in dXi:
        knots.append(knots[-1] + dx)

    # cumulative knots = [0,1,3,6,7]

    inner = knots[1:-1]

    Xi = [0.0]*(p+1) + inner + [knots[-1]]*(p+1)

    return Xi, knots


Xi, knots = build_knot_vector(dXi, p)

# Xi = [0,0,0,0,1,3,6,7,7,7,7]
# knots = [0,1,3,6,7]


# =====================================================
# Cox-de Boor recursion
# =====================================================
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
        left = ((x - Xi[i]) / d1) * N(i, p-1, x, Xi)

    if d2 > 0:
        right = ((Xi[i+p+1] - x) / d2) * N(i+1, p-1, x, Xi)

    return left + right


# =====================================================
# Number of basis functions
# =====================================================
n_basis = len(Xi) - p - 1


# =====================================================
# Write data
# =====================================================
xmin = Xi[0]
xmax = Xi[-1]

with open("bspline.dat", "w") as f:

    for x in np.linspace(xmin, xmax, 1000):

        vals = [x]

        for i in range(n_basis):
            vals.append(N(i, p, x, Xi))

        f.write(" ".join(map(str, vals)) + "\n")


# =====================================================
# X-axis labels only at knot locations
# (excluding 0 as requested)
# =====================================================
tick_values = knots[1:]   # [1,3,6,7]

tick_string = ", ".join([f'"{v}" {v}' for v in tick_values])

caption = "Polynomial order: p = 3, knot spans: ΔΞ = {1,2,3,1}"


# =====================================================
# GNUplot
# =====================================================
gp = f"""
set terminal qt size 1000,650
set title "Cubic B-Spline Basis Functions"
set xlabel "Ξ"
set ylabel "N"
set key outside
unset grid

set xrange [{xmin}:{xmax}]
set yrange [0:1.05]

# only selected ticks
set xtics ({tick_string})

# vertical dotted lines
set arrow from 1,0 to 1,1 nohead dt 2 lw 1
set arrow from 3,0 to 3,1 nohead dt 2 lw 1
set arrow from 6,0 to 6,1 nohead dt 2 lw 1
set arrow from 7,0 to 7,1 nohead dt 2 lw 1

set label 1 "{caption}" at screen 0.5,0.03 center

plot for [i=2:{n_basis+1}] 'bspline.dat' using 1:i with lines lw 2 title sprintf("N_%d,3", i-2)

pause -1
"""

subprocess.run(["gnuplot", "-persist"], input=gp, text=True)