#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
import numpy as np
from math import comb


#  PART 1 (background/reuse): Cox-de Boor + element Bezier
#  extraction, only used here to *generate* a realistic test case.

def bspline_basis(i, p, U, xi):
    if p == 0:
        if (U[i] <= xi < U[i+1]) or (xi == U[-1] and U[i] <= xi <= U[i+1]):
            return 1.0
        return 0.0
    left = 0.0
    right = 0.0
    if U[i+p] != U[i]:
        left = ((xi - U[i]) / (U[i+p] - U[i])) * bspline_basis(i, p-1, U, xi)
    if U[i+p+1] != U[i+1]:
        right = ((U[i+p+1] - xi) / (U[i+p+1] - U[i+1])) * bspline_basis(i+1, p-1, U, xi)
    return left + right

def bernstein(p, xi):
    return np.array([comb(p, i) * (xi**i) * ((1-xi)**(p-i)) for i in range(p+1)])

def element_extraction_operator(U, p, element):
    xi_left = U[element+p]
    xi_right = U[element+p+1]
    sample = np.linspace(0.1, 0.9, p+1)
    B, N = [], []
    for s in sample:
        x = xi_left + s*(xi_right-xi_left)
        B.append(bernstein(p, s))
        N.append([bspline_basis(i, p, U, x) for i in range(element, element+p+1)])
    B, N = np.array(B), np.array(N)
    return np.linalg.solve(B, N).T          # C[local_i, bernstein_a]

def bezier_subdivide(C, alpha):
    C = list(np.asarray(C, dtype=float))
    p = len(C) - 1
    left = [C[0]]
    right = [C[-1]]
    pts = C[:]
    for _ in range(p):
        pts = [(1-alpha)*pts[k] + alpha*pts[k+1] for k in range(len(pts)-1)]
        left.append(pts[0])
        right.append(pts[-1])
    right.reverse()
    return np.array(left), np.array(right)



def diff_at_1(C, p, d):
    return sum((-1)**j * comb(d, j) * C[p - j] for j in range(d + 1))

def diff_at_0(C, p, d):
    return sum((-1)**(d - j) * comb(d, j) * C[j] for j in range(d + 1))


def continuity_and_ratio(CA, CB, p, tol=1e-8):
    CA = np.asarray(CA, dtype=float)
    CB = np.asarray(CB, dtype=float)

    # d = 0 must already match (shared value at the boundary)
    d0_err = abs(diff_at_1(CA, p, 0) - diff_at_0(CB, p, 0))
    if d0_err > 1e-6:
        raise ValueError(f"Not C0 continuous at this boundary")

    DA1 = diff_at_1(CA, p, 1)
    DB1 = diff_at_0(CB, p, 1)

    if abs(DB1) < tol and abs(DA1) < tol:
        r = 1.0   # degenerate case (both first differences vanish)
    else:
        r = DB1 / DA1

    dmax = 1
    for d in range(2, p + 1):
        DAd = diff_at_1(CA, p, d)
        DBd = diff_at_0(CB, p, d)
        lhs = DBd
        rhs = (r ** d) * DAd
        if abs(lhs - rhs) < tol * max(1.0, abs(lhs), abs(rhs)):
            dmax = d
        else:
            break

    multiplicity = p - dmax
    alpha = 1.0 / (r + 1.0)
    return dmax, multiplicity, r, alpha


def bezier_rows_to_knot_spans(C_list, p, tol=1e-8, verbose=True):
    n = len(C_list)
    if n < p + 1:
        raise ValueError(f"Need at least p+1={p+1} row operators, got {n}")

    r_list, mult_list, dmax_list = [], [], []
    for i in range(n - 1):
        dmax, mult, r, alpha = continuity_and_ratio(C_list[i], C_list[i+1], p, tol)
        r_list.append(r)
        mult_list.append(mult)
        dmax_list.append(dmax)
        if verbose:
            tag = ("ZERO knot-multiplicity -> artificial split (merge)"
                   if mult == 0 else f"real knot, multiplicity = {mult}")
            print(f"  boundary {i}-{i+1}:  dmax={dmax},  r=s{i+1}/s{i}={r: .6f}   [{tag}]")

    # relative element sizes s_0..s_{n-1} (absolute scale is not
    # recoverable from extraction operators alone -> fix s_0 = 1)
    s = [1.0]
    for r in r_list:
        s.append(s[-1] * r)
    s = np.array(s)

    # merge across artificial (zero-multiplicity) boundaries
    merged = [s[0]]
    for i, mult in enumerate(mult_list):
        if mult == 0:
            merged[-1] += s[i + 1]
        else:
            merged.append(s[i + 1])
    merged = np.array(merged)

    if len(merged) != p + 1:
        print(f"  WARNING: expected {p+1} real knot spans, found {len(merged)}.")

    normalized_spans = merged / merged.sum()
    return normalized_spans, mult_list, r_list




if __name__ == "__main__":
    p = 2
    U = [0, 0, 0, 1, 2, 3, 3, 3]          # spans [0,1],[1,2],[2,3] -> equal spans
    num_elements = len(U) - 1 - 2*p

    g = 2
    true_rows = []
    for e in range(num_elements):
        i_local = g - e
        Ce = element_extraction_operator(U, p, e)
        true_rows.append(Ce[i_local, :])

    print("True (un-split) row operators for basis function g =", g)
    for e, row in enumerate(true_rows):
        print(f"  element {e}: {np.round(row, 6)}")

    print("\nTrue knot spans (normalized): [1/3, 1/3, 1/3] "
          "since U spans are [0,1],[1,2],[2,3] (all length 1)\n")

    left_half, right_half = bezier_subdivide(true_rows[1], alpha=0.5)

    split_rows = [true_rows[0], left_half, right_half, true_rows[2]]

    print("Row operators AFTER an artificial split of element 1 (n=4 rows):")
    for e, row in enumerate(split_rows):
        print(f"  row {e}: {np.round(row, 6)}")

    print("\nRunning bezier_rows_to_knot_spans on the SPLIT data "
          "(should detect the fake boundary and recover 3 equal spans):\n")
    spans, mult_list, r_list = bezier_rows_to_knot_spans(split_rows, p)

    print("\nRecovered normalized knot spans:", np.round(spans, 6))
    print("Detected multiplicities at each of the", len(mult_list), "boundaries:", mult_list)
