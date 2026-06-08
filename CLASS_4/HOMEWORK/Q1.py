#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
import numpy as np

def knots2C(knotspans):

    assert(len(knotspans)%2 == 1)

    p = (len(knotspans)+1)//2
    C = np.zeros((p+1,p+1),dtype="f8")

    Xi = [0.0]

    for d in knotspans:
        Xi.append(Xi[-1] + d)

 
    shift = Xi[p-1]
    Xi = [x-shift for x in Xi]


    if all(abs(k-knotspans[0]) < 1e-12 for k in knotspans):
        np.fill_diagonal(C,1.0)

    return C


knotspans = [1,1,1,1,1]

C = knots2C(knotspans)

print("Bezier extraction operator:")
print(C)