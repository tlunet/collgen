#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 19:05:48 2022

@author: telu
"""
import os
import sys
import sympy as sy

sys.path.append(os.path.realpath('../'))
from collgen.thomas import genCollFromRK


nodes, Q, weights = genCollFromRK('RK4')


dt, f = sy.symbols(r'\Delta{t} f')

coll = sy.eye(5) - dt*Q*f

coll[-1, 0] = -sy.Rational(1, 6)*dt*f
coll[-1, -2] = -sy.Rational(1, 6)*dt*f
coll[1, 0] = -sy.Rational(1, 2)*dt*f
coll[2, 1] = -sy.Rational(1, 2)*dt*f
coll[3, 2] = -dt*f
coll[-1, 1] = -sy.Rational(1, 3)*dt*f
coll[-1, 2] = -sy.Rational(1, 3)*dt*f

# Zero-to-node formulation of a collocation problem !

sy.print_latex(coll[:-1, :-1])
