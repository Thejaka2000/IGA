#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
import numpy as np
from math import comb

def splitOperators(p, alpha):

    assert(0 <= alpha)
    assert(alpha <= 1)

    SA = np.zeros((p+1, p+1), dtype="f8")
    SB = np.zeros((p+1, p+1), dtype="f8")

    # Left split operator
    for a in range(p+1):
        for b in range(a, p+1):

            SA[a,b] = (
                comb(b, a)
                * alpha**a
                * (1-alpha)**(b-a)
            )

    # Right split operator
    for a in range(p+1):
        for b in range(a+1):

            SB[a,b] = (
                comb(p-b, a-b)
                * alpha**(a-b)
                * (1-alpha)**(p-a)
            )

    return SA, SB

p = 3
alpha = 0.5

SA, SB = splitOperators(p, alpha)

print("SA =")
print(SA)

print("\nSB =")
print(SB)