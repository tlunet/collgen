#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description
-----------

Module containing a function generating Collocation matrix
following Thomas'implementation of the RungeKutta sweeper in pySDC.
"""
import numpy as np

from .butcher import tablesExplicit, tablesImplicit

def genCollFromRK(scheme):
    # Get Butcher coefficients
    if scheme in tablesExplicit:
        table = tablesExplicit[scheme]
    elif scheme in tablesImplicit:
        table = tablesImplicit[scheme]
    else:
        raise ValueError(f'RK scheme {scheme} not implemented')

    A = np.array(table['A'])
    b = np.array(table['b'])
    c = np.array(table['c'])

    # Q-matrix formulation
    M = b.size + 1
    nodes = np.concatenate((c, [1.]))
    Q = np.zeros((M, M))
    Q[:-1, :-1] = A
    Q[-1, :-1] = b
    weights = b

    return nodes, Q, weights