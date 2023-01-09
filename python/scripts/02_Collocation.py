#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 12:56:34 2023

@author: cpf5546
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.realpath('../'))

from collgen.qmatrix import genCollocation
from collgen.vectorize import matVecMul
from collgen.plots import plotAccuracyContour

M = 4
distr = 'LEGENDRE'
quadType = 'LOBATTO'
nodes, Q, weights = genCollocation(M, distr, quadType)

# Picard iteration
reLam = np.linspace(-3, 0.5, 501)
imLam = np.linspace(-3, 3, 500)
lam = reLam[:, None] + 1j*imLam[None, :]
lamDt = np.ravel(lam)[None, :]

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
plt.savefig(f'Collocation_M{M}_{quadType}_{distr}_Picard.svg')
