#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 19:05:48 2022

@author: telu
"""
import os
import sys

import numpy as np
import sympy as sy

sys.path.append(os.path.realpath('../'))
from collgen.thomas import genCollFromRK
from collgen.vectorize import matVecMul
from collgen.plots import plotAccuracyContour

nodes, Q, weights = genCollFromRK('RK4')

# Some sympy representation
dt, f = sy.symbols(r'\Delta{t} f')

coll = sy.eye(5) - dt*Q*f

coll[-1, 0] = -sy.Rational(1, 6)*dt*f
coll[-1, -2] = -sy.Rational(1, 6)*dt*f
coll[1, 0] = -sy.Rational(1, 2)*dt*f
coll[2, 1] = -sy.Rational(1, 2)*dt*f
coll[3, 2] = -dt*f
coll[-1, 1] = -sy.Rational(1, 3)*dt*f
coll[-1, 2] = -sy.Rational(1, 3)*dt*f

# sy.print_latex(coll[:-1, :-1])

# Picard iteration
reLam = np.linspace(-4, 0.5, 501)
imLam = np.linspace(-3, 3, 500)
lam = reLam[:, None] + 1j*imLam[None, :]
lamDt = np.ravel(lam)[None, :]

# Q = Q[:-1, :-1]

Qvec = lamDt*Q[..., None]
Qvec = Qvec.transpose((2,0,1))

nIter = 5
u0 = np.ones(Q.shape[0])

uPic = u0
for k in range(nIter):
    uPic = u0 + matVecMul(Qvec, uPic)

uExact = np.exp(lamDt).ravel()

err = np.abs(uExact-uPic[:, -1])  # Comparison with last node

stab = np.abs(uPic[:, -1]).reshape(lam.shape)
errMax = err.reshape(lam.shape)

# Plot discretization error on complex plane
plotAccuracyContour(reLam, imLam, errMax, stab)
