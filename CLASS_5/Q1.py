#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Created: by Thejaka Jayasinghe
import numpy

def recoverKnotSpans(CRs):
    assert(len(CRs) > 0)

    p = len(CRs[0]) - 1

    knot_spans = numpy.empty((p + 1, p + 1), dtype="f8")

    # First knot span [0,1]
    knot_spans[:, 0] = 0.0
    knot_spans[:, -1] = 1.0

    for j in range(1, p):
        knot_spans[:, j] = CRs[j - 1]

    return knot_spans

CRs = [
    numpy.array([1.0, 0.5, 0.0]),
    numpy.array([1.0, 0.25, 0.0])
]

knot_spans = recoverKnotSpans(CRs)

print(knot_spans)