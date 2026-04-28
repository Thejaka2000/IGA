#!/usr/bin/env python3
import numpy as np
import subprocess

#  INPUT
dXi = [1,1,1]
p = 2   # quadratic

# Build knot vector
def build_knot_vector(dXi, p):

    knots = [0.0]

    for dx in dXi:
        knots.append(knots[-1] + dx)

    inner = knots[1:-1]

    Xi = [0.0]*(p+1) + inner + [knots[-1]]*(p+1)

    return Xi

Xi = build_knot_vector(dXi,p)

# Cox-de Boor recursion
def N(i,p,x,Xi):

    if p == 0:
        if (Xi[i] <= x < Xi[i+1]) or (x == Xi[-1] and Xi[i] <= x <= Xi[i+1]):
            return 1.0
        return 0.0

    left = 0.0
    right = 0.0

    d1 = Xi[i+p] - Xi[i]
    d2 = Xi[i+p+1] - Xi[i+1]

    if d1 > 0:
        left = ((x-Xi[i])/d1)*N(i,p-1,x,Xi)

    if d2 > 0:
        right = ((Xi[i+p+1]-x)/d2)*N(i+1,p-1,x,Xi)

    return left + right

# Number of basis functions
n_basis = len(Xi)-p-1

# quadratic case with 5 basis functions
# pair ends, keep middle

mid = n_basis // 2

xmin = Xi[0]
xmax = Xi[-1]

with open("support.dat","w") as f:

    for x in np.linspace(xmin,xmax,600):

        vals = [x]

        # left pair
        vals.append(N(0,p,x,Xi) + N(1,p,x,Xi))

        # middle
        vals.append(N(mid,p,x,Xi))

        # right pair
        vals.append(N(mid+1,p,x,Xi) + N(mid+2,p,x,Xi))

        f.write(" ".join(map(str,vals)) + "\n")

# GNUplot
gp = f"""
set terminal qt size 900,600
set title "Supporting Functions"
set xlabel "Ξ"
set ylabel "N"
set grid
set key outside
set xrange [{xmin}:{xmax}]
set yrange [0:1.05]

plot \
'support.dat' using 1:2 w l lw 3 title 'N0+N1', \
'support.dat' using 1:3 w l lw 3 title 'N2', \
'support.dat' using 1:4 w l lw 3 title 'N3+N4'

pause -1
"""

subprocess.run(["gnuplot","-persist"], input=gp, text=True)
