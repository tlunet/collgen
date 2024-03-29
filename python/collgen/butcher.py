#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description
-----------

Module containing butcher tables for Runge-Kutta methods

References
----------
.. [1] Wang, R., & Spiteri, R. J. (2007). Linear instability of the fifth-order
       WENO method. SIAM Journal on Numerical Analysis, 45(5), 1871-1901.
.. [2] Alexander, R. (1977). Diagonally implicit Runge-Kutta methods for stiff
       ODE's. SIAM Journal on Numerical Analysis, 14(6), 1006-1021.
.. [3] Wanner, G., & Hairer, E. (1996). Solving ordinary differential equations
       II. Springer Berlin Heidelberg.
.. [4] Butcher, J.C. (2003). Numerical methods for Ordianry Differential
       Equations. John Wiley & Sons.
"""

tablesExplicit = \
    {'RK4':  # Classical Runge Kutta method of order 4 (cf Wikipedia)
     {'A': [[0, 0, 0, 0],
            [0.5, 0, 0, 0],
            [0, 0.5, 0, 0],
            [0, 0, 1, 0]],
      'b': [1/6, 1/3, 1/3, 1/6],
      'c': [0, 1/2, 1/2, 1],
      'c1': -1/120, 'lowStorage': True, 'order': 4},
     'RK4-2':  # The 3/8-rule due to Kutta, order 4  (cf Wikipedia)
     {'A': [[0, 0, 0, 0],
            [1/3, 0, 0, 0],
            [-1/3, 1, 0, 0],
            [1, -1, 1, 0]],
      'b': [1/8, 3/8, 3/8, 1/8],
      'c': [0, 1/3, 2/3, 1],
      'c1': -1/120, 'lowStorage': True, 'order': 4},
     'FE':  # Forward Euler method (cf Wikipedia)
     {'A': [[0]],
      'b': [1],
      'c': [0],
      'c1': -1/2, 'lowStorage': True, 'order': 1},
     'RK53':  # Method in 5 steps of order 3 from Wang & Spiteri [1]
     {'A': [[0, 0, 0, 0, 0],
            [1/7, 0, 0, 0, 0],
            [0, 3/13, 0, 0, 0],
            [0, 0, 1/3, 0, 0],
            [0, 0, 0, 2/3, 0]],
      'b': [1/4, 0, 0, 0, 3./4],
      'c': [0, 1/7, 3/16, 1/3, 2/3],
      'c1': -1/312, 'lowStorage': True, 'order': 3},
     'RK21':  # Method in 2 steps of order 1 from Wang & Spiteri [1]
     {'A': [[0, 0],
            [3./4, 0]],
      'b': [-1./3, 4./3],
      'c': [0, 3./4],
      'c1': 1/2, 'lowStorage': True, 'order': 1},
     'RK2':  # Classical Runge Kutta method of order 2 (cf Wikipedia)
     {'A': [[0, 0],
            [1./2, 0]],
      'b': [0, 1.],
      'c': [0, 1./2],
      'c1': -1/6, 'lowStorage': True, 'order': 2},
     'HEUN2':  # Heun method of order 2 (cf Wikipedia)
     {'A': [[0, 0],
            [1., 0]],
      'b': [1./2, 1./2],
      'c': [0, 1.],
      'c1': -1/6, 'lowStorage': True, 'order': 2},
     'RK32':  # Method in 3 steps of order 2 from Wang & Spiteri [1]
     {'A': [[0, 0, 0],
            [1./3, 0, 0],
            [0, 1., 0]],
      'b': [1./2, 0, 1./2],
      'c': [0, 1./3, 1.],
      'c1': -1/24, 'lowStorage': True, 'order': 2},
     'RK33':  # Method in 3 steps of order 3 from Wang & Spiteri [1]
     {'A': [[0, 0, 0],
            [1., 0, 0],
            [1./4, 1./4, 0]],
      'b': [1./6, 1./6, 2./3],
      'c': [0, 1., 1./2],
      'c1': -1/24, 'lowStorage': False, 'order': 3},
     'RK65':  # Method in 6 steps of order 5, (236a) from Butcher [4]
     {'A': [[0, 0, 0, 0, 0, 0],
            [0.25, 0, 0, 0, 0, 0],
            [1./8, 1./8, 0, 0, 0, 0],
            [0, 0, 0.5, 0, 0, 0],
            [3./16, -3./8, 3./8, 9./16, 0, 0],
            [-3./7, 8./7, 6./7, -12./7, 8./7, 0]],
      'b': [7./90, 0, 32./90, 12./90, 32./90, 7./90],
      'c': [0, 0.25, 0.25, 0.5, 0.75, 1.],
      'c1': -7/11520, 'lowStorage': False, 'order': 5},
     }

tablesImplicit = \
    {'BE':  # Backward Euler method (also SDIRK1, see [2])
     {'A': [[1.]],
      'b': [1.],
      'c': [1.],
      'c1': 1/2, 'order': 1},
     'TRAPZ':  # Trapeze method (cf Wikipedia)
     {'A': [[0, 0],
            [1./2, 1./2]],
      'b': [1./2, 1./2],
      'c': [0., 1.],
      'c1': 1/12, 'order': 2},
     'GAUSS-LG':  # Gauss-Legendre method of order 4 (cf Wikipedia)
     {'A': [[0.25, 0.25-1./6*3**(0.5)],
            [0.25+1./6*3**(0.5), 0.25]],
      'b': [0.5, 0.5],
      'c': [0.5-1./6*3**(0.5), 0.5+1./6*3**(0.5)],
      'c1': -1/720, 'order': 4},
     'SDIRK2':  # First S-stable Diagonally Implicit Runge Kutta method of
                # order 2 in two stages, from Alexander [2]
     {'A': [[1-2**0.5/2, 0],
            [2**0.5/2, 1-2**0.5/2]],
      'b': [2**0.5/2, 1-2**0.5/2],
      'c': [1-2**0.5/2, 1.],
      'c1': 0.0404401145198803, 'order': 2},
     'SDIRK2-2':  # Second S-stable Diagonally Implicit Runge Kutta method of
                  # order 2 in two stages, from Alexander [2]
     {'A': [[1+2**0.5/2, 0],
            [-2**0.5/2, 1+2**0.5/2]],
      'b': [-2**0.5/2, 1+2**0.5/2],
      'c': [1+2**0.5/2, 1.],
      'c1': -1.37377344785315, 'order': 2},
     'SDIRK3':  # S-stable Diagonally Implicit Runge Kutta method of
                # order 3 in three stages, from Alexander [2]
     {'A': [[0.43586652150845967, 0, 0],
            [0.28206673924577014, 0.43586652150845967, 0],
            [1.2084966491760119, -0.6443631706844715, 0.43586652150845967]],
      'b': [1.2084966491760119, -0.6443631706844715, 0.43586652150845967],
      'c': [0.43586652150845967, 0.7179332607542298, 1.],
      'c1': -0.0258970846506337, 'order': 3},
     'SDIRK54':  # S-stable Diagonally Implicit Runge Kutta method of
                 # order 4 in five stages, from Wanner and Hairer [3]
     {'A': [[1/4, 0, 0, 0, 0],
            [1/2, 1/4, 0, 0, 0],
            [17/50, -1/25, 1/4, 0, 0],
            [371/1360, -137/2720, 15/544, 1/4, 0],
            [25/24, -49/48, 125/16, -85/12, 1/4]],
      'b': [25/24, -49/48, 125/16, -85/12, 1/4],
      'c': [1/4, 3/4, 11/20, 1/2, 1],
      'c1': -8.46354080300333e-4, 'order': 4},
     }
