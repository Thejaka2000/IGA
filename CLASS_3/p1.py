#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
# $Author$
# $Date$
# $URL$
__giturl__ = "$URL$"


from itertools import product


def generate_derivative_table(npd, d):

    combinations = []

    def compositions(remaining, current):
        if len(current) == npd - 1:
            combinations.append(current + [remaining])
            return

        for value in range(remaining, -1, -1):
            compositions(remaining - value, current + [value])

    compositions(d, [])

    table = []

    for deriv in combinations:

        bits = []

        for i, count in enumerate(deriv):

            bits.extend("0" * count)

            if i != npd - 1:
                bits.append("1")

        binary = "".join(bits)

        decimal = int(binary, 2)

        table.append((tuple(deriv), binary, decimal))

    return table


def print_table(npd, d):

    print(f"\nnpd = {npd}, d = {d}")
    print("-" * 40)
    print(f"{'Derivative':15s}{'Binary':10s}{'s'}")

    for deriv, binary, decimal in generate_derivative_table(npd, d):
        print(f"{str(deriv):15s}{binary:10s}{decimal}")


# Ouput for each case
print_table(1, 2)
print_table(2, 1)
print_table(2, 2)
print_table(3, 2)