#!/usr/bin/env python3
# Multinomial coefficients example for p = 2 , n = 2
# Expansion of:
# (x0 + x1 + x2)^2
#
# Parametric dimension n = 2 means variables:
# x0, x1, x2

import math
from itertools import product
import sympy as sp

# ---------------------------------
# symbols
# ---------------------------------
x0, x1, x2 = sp.symbols('x0 x1 x2')

p = 2
n = 2

vars_list = [x0, x1, x2]

# ---------------------------------
# Find all tuples (m0,m1,m2)
# such that m0+m1+m2 = p
# ---------------------------------
terms = []

for m in product(range(p + 1), repeat=n + 1):
    if sum(m) == p:
        terms.append(m)

# ---------------------------------
# Print basis terms
# ---------------------------------
print(f"Expansion of (x0 + x1 + x2)^{p}\n")

expr = 0

for m in terms:
    m0, m1, m2 = m

    coeff = math.factorial(p) / (
        math.factorial(m0) *
        math.factorial(m1) *
        math.factorial(m2)
    )

    basis = coeff * x0**m0 * x1**m1 * x2**m2

    print(f"m = {m}   coefficient = {int(coeff)}   term = {basis}")

    expr += basis

# ---------------------------------
# Final expression
# ---------------------------------
print("\nExpanded polynomial:\n")
print(sp.expand(expr))