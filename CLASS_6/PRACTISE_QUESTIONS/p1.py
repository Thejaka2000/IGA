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


# Bernstein -> Power basis transformation matrix


def bernstein_to_power_matrix(p):

    M = np.zeros((p + 1, p + 1), dtype=int)

    # Bernstein polynomial B_i^p
    for i in range(p + 1):

        # Expand (1-t)^(p-i)
        for k in range(p - i + 1):

            coeff = comb(p, i) * comb(p - i, k) * (-1) ** k

            power = i + k

            M[i, power] = coeff

    return M



# Print matrices


for p in range(1, 6):

    print("=" * 50)
    print(f"Degree p = {p}\n")

    M = bernstein_to_power_matrix(p)

    print(M)
    print()
