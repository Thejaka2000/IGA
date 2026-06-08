#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
import numpy as np
from math import comb

# Split Operator
def splitOperators(p, alpha):

    SA = np.zeros((p+1,p+1))
    SB = np.zeros((p+1,p+1))

    for a in range(p+1):
        for b in range(a,p+1):

            SA[a,b] = (
                comb(b,a)
                * alpha**a
                * (1-alpha)**(b-a)
            )

    for a in range(p+1):
        for b in range(a+1):

            SB[a,b] = (
                comb(p-b,a-b)
                * alpha**(a-b)
                * (1-alpha)**(p-a)
            )

    return SA,SB

# Cubic Bezier
p = 3

P = np.array([
    [0.0,0.0],
    [1.0,2.0],
    [3.0,2.0],
    [4.0,0.0]
])

alpha = 0.5
SA,SB = splitOperators(p,alpha)

# New control points
PA = SA.T @ P
PB = SB.T @ P

# Bernstein basis
def bernstein(p,i,t):
    return comb(p,i)*(t**i)*(1-t)**(p-i)

def bezier(P,t):

    x = 0.0
    y = 0.0

    for i in range(len(P)):

        B = bernstein(len(P)-1,i,t)

        x += B*P[i,0]
        y += B*P[i,1]

    return x,y

# Original curve
with open("original.dat","w") as f:
    for t in np.linspace(0,1,400):
        x,y = bezier(P,t)
        f.write(f"{x} {y}\n")

# Refined curve
with open("refined.dat","w") as f:
    for t in np.linspace(0,1,200):
        x,y = bezier(PA,t)
        f.write(f"{x} {y}\n")
    for t in np.linspace(0,1,200):
        x,y = bezier(PB,t)
        f.write(f"{x} {y}\n")
print("Files written")