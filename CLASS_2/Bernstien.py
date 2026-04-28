#!/usr/bin/env python3

import math

def cal_B(p, t):
    if not (0 <= t <= 1):
        raise ValueError("t must satisfy 0 <= t <= 1")

    B = []

    for i in range(p + 1):
        val = math.comb(p, i) * (t ** i) * ((1 - t) ** (p - i))
        B.append(val)

    return B


# Example
p = int(input("Enter order p: "))
t = float(input("Enter value t (0 <= t <= 1): "))

B = cal_B(p, t)

print("p =", p)
print("t =", t)
print("Bernstein basis values =", B)
print("Sum =", sum(B))
