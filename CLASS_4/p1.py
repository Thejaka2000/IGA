#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
# $Author$
# $Date$
# $URL$
__giturl__ = "$URL$"

import numpy as np
from math import comb


def bspline_basis(i, p, U, xi):
    if p == 0:
        if (U[i] <= xi < U[i+1]) or (xi == U[-1] and U[i] <= xi <= U[i+1]):
            return 1.0
        return 0.0
    left = 0.0
    right = 0.0
    if U[i+p] != U[i]:
        left = ((xi - U[i]) / (U[i+p] - U[i])) * \
               bspline_basis(i, p-1, U, xi)
    if U[i+p+1] != U[i+1]:
        right = ((U[i+p+1] - xi) / (U[i+p+1] - U[i+1])) * \
                bspline_basis(i+1, p-1, U, xi)
    return left + right


# General Bernstein basis (works for ANY degree p)
# B_{i,p}(xi) = C(p,i) * xi^i * (1-xi)^(p-i),  i = 0..p
def bernstein(p, xi):
    return np.array([
        comb(p, i) * (xi**i) * ((1 - xi)**(p - i))
        for i in range(p + 1)
    ])


# Compute extraction operator for ONE element

def extraction_operator(U, p, element):
    xi_left = U[element + p]
    xi_right = U[element + p + 1]
    sample = np.linspace(0.1, 0.9, p + 1)  # p+1 sample points -> square system

    B = []
    N = []
    for s in sample:
        # map local coord [0,1] -> physical knot span
        x = xi_left + s * (xi_right - xi_left)
        B.append(bernstein(p, s))
        row = [bspline_basis(i, p, U, x) for i in range(element, element + p + 1)]
        N.append(row)

    B = np.array(B)
    N = np.array(N)
    # Solve N = B C^T  ->  C^T = B^-1 N  ->  C = (B^-1 N)^T
    C = np.linalg.solve(B, N).T
    return C


# Extraction operator of every element
def run_case(p, interior_knots, domain_end):
    # open knot vector: (p+1) zeros, strictly-interior knots, (p+1) copies of domain_end
    U = [0]*(p+1) + list(interior_knots) + [domain_end]*(p+1)
    U = np.array(U, dtype=float)

    num_elements = len(U) - 1 - 2*p
    print("="*60)
    print(f"Degree p = {p}")
    print("Knot vector U =", U.tolist())
    print("Number of elements =", num_elements)

    operators = []
    for e in range(num_elements):
        C = extraction_operator(U, p, e)
        operators.append(C)
        print(f"\nElement {e}  (span [{U[e+p]}, {U[e+p+1]}])")
        print(np.round(C, 6))
    return U, operators


# Run for degree 2, 3, 4

if __name__ == "__main__":
    # degree 2: U = [0,0,0,1,2,3,3,3]  -> spans [0,1],[1,2],[2,3]
    run_case(p=2, interior_knots=[1, 2], domain_end=3)

    # degree 3: U = [0,0,0,0,1,2,3,3,3,3]
    run_case(p=3, interior_knots=[1, 2], domain_end=3)

    # degree 4: U = [0,0,0,0,0,1,2,3,3,3,3,3]
    run_case(p=4, interior_knots=[1, 2], domain_end=3)
