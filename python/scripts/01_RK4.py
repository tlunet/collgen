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
