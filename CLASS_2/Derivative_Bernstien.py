#!/usr/bin/env python3
import math
import subprocess

# Bernstein basis and derivatives at any t
# B_p,i^(d)(t)
# d = derivative order
def bernstein_derivative(p, i, d, t):
    """
    Compute d-th derivative of Bernstein basis B_i^p(t)
    """

    if d == 0:
        return math.comb(p, i) * (t**i) * ((1 - t)**(p - i))

    if d > p:
        return 0.0

    coeff = math.factorial(p) / math.factorial(p - d)

    val = 0.0

    for k in range(d + 1):
        j = i - k

        if 0 <= j <= (p - d):
            val += ((-1)**(d - k)) * math.comb(d, k) * \
                   math.comb(p - d, j) * \
                   (t**j) * ((1 - t)**((p - d) - j))

    return coeff * val


# Get all basis values for given derivative d
def cal_B(p, d, t):

    if not (0 <= t <= 1):
        raise ValueError("t must satisfy 0 <= t <= 1")

    B = []

    for i in range(p + 1):
        B.append(bernstein_derivative(p, i, d, t))

    return B


# Print endpoint values
def print_endpoints(p, max_d):

    print(f"\nBernstein basis values for p = {p}\n")

    for d in range(max_d + 1):

        print(f"Derivative order d = {d}")

        B0 = cal_B(p, d, 0.0)
        B1 = cal_B(p, d, 1.0)

        print("t = 0 :", B0)
        print("t = 1 :", B1)
        print()


# Plot basis curves using gnuplot
def plot_basis(p):

    fname = "bernstein.dat"

    with open(fname, "w") as f:

        npts = 200

        for k in range(npts + 1):

            t = k / npts
            vals = cal_B(p, 0, t)

            line = [t] + vals
            f.write(" ".join(map(str, line)) + "\n")

    gp = f"""
    set title 'Bernstein Basis Functions p={p}'
    set xlabel 't'
    set ylabel 'B_i^p(t)'
    set grid
    plot for [i=2:{p+2}] '{fname}' using 1:i with lines lw 2 title sprintf('B_%d',i-2)
    pause -1
    """

    subprocess.run(["gnuplot", "-persist"], input=gp, text=True)


# -------------------------------------------------
# MAIN
# -------------------------------------------------
p = int(input("Enter polynomial order p = "))
max_d = int(input("Enter max derivative order d = "))

print_endpoints(p, max_d)

plot_basis(p)
