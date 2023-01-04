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
import matplotlib.pyplot as plt

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
reLam = np.linspace(-3, 0.5, 501)
imLam = np.linspace(-3, 3, 500)
lam = reLam[:, None] + 1j*imLam[None, :]
lamDt = np.ravel(lam)[None, :]

# Q = Q[:-1, :-1]

Qvec = lamDt*Q[..., None]
Qvec = Qvec.transpose((2,0,1))

nIter = 4
u0 = np.ones(Q.shape[0])
uPic = np.zeros((nIter+1, Qvec.shape[0], u0.size), dtype=Qvec.dtype)

uPic[0] = u0
for k in range(nIter):
    uPic[k+1] = u0 + matVecMul(Qvec, uPic[k])

uExact = np.exp(lamDt).ravel()

err = np.abs(uExact-uPic[:, :, -1])  # Comparison with last node

stab = np.abs(uPic[:, :, -1]).reshape((nIter+1, *lam.shape))
errMax = err.reshape((nIter+1, *lam.shape))

# Plot discretization error on complex plane
fig, axes = plt.subplots(1, 4)
axes = axes.ravel()
for k in range(nIter):
    plotAccuracyContour(reLam, imLam, errMax[k+1], stab[k+1], axe=axes[k])
fig.set_size_inches(12.29,  3.15)
fig.tight_layout()
fig.tight_layout()
fig.tight_layout()
plt.savefig('RK4_Picard.svg')
