#!/usr/bin/env python3
import math
import sympy as sp

# Symbolic parameter
xi = sp.Symbol('xi')

# -----------------------------
# Bernstein basis function
# -----------------------------
def bernstein_basis_1D(p):
    basis = []

    for i in range(p + 1):
        coeff = math.factorial(p) / (
            math.factorial(i) * math.factorial(p - i)
        )

        B = coeff * (1 - xi)**(p - i) * xi**i
        basis.append(sp.expand(B))

    return basis


# -----------------------------
# Main
# -----------------------------
p = int(input("Enter polynomial order p = "))

basis = bernstein_basis_1D(p)

print(f"\n1D Bernstein Basis Functions of order {p}:\n")

for i, B in enumerate(basis):
    print(f"B{i}^{p}(xi) = {B}")