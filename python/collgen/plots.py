#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 12:13:43 2023

@author: telu
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

def plotAccuracyContour(reLam, imLam, err, stab=None,
                        eMin=-7, eMax=0, nLevels=22,
                        figName=None, axe=None):
    """
    2D contour plot of an error given in parameters for many complex values

    Parameters
    ----------
    reLam : 1darray (nR,)
        The values for real part of lambda
    imLam : 1darray (nI,)
        The values for imaginary part of lambda.
    err : 2darray (nR, nI)
        The error values for each lambda
    stab : 2darray (nR, nI), optional
        Amplification factor associated to the error. The default is None.
    eMin : int, optional
        Minimum exponent to be shown for the error. The default is -7.
    eMax : int, optional
        Maximum exponent to be shown for the error. The default is 0.
    nLevels : int, optional
        Number of level to show on the contour plot. The default is 22.
    figName : str, optional
        Name for the generated figure (ignored if axe is given).
        The default is None.
    axe : plt.AxesSubplot, optional
        The axes to plot the figure on (by default, uses plt.gca())
    """
    coords = np.meshgrid(reLam.ravel(), imLam.ravel(), indexing='ij')
    levels = np.logspace(eMin, eMax, num=nLevels)
    err[err < 10**eMin] = 10**eMin
    err[err > 10**eMax] = 10**eMax
    ticks = [10**(i) for i in range(eMin, eMax+1)]

    if axe is None:
        fig = plt.figure(figName)
        axe = plt.gca()
    else:
        fig = axe.get_figure()

    cm = axe.contourf(*coords, err, levels=levels, locator=ticker.LogLocator())
    fig.colorbar(cm, ticks=ticks, format=ticker.LogFormatter())
    axe.contour(*coords, err, levels=ticks,
                colors='k', linestyles='--', linewidths=0.75)
    axe.contour(*coords, stab, levels=[1], colors='gray')
    axe.hlines(0, coords[0].min(), coords[0].max(),
               colors='black', linestyles='--')
    axe.vlines(0, coords[1].min(), coords[1].max(),
               colors='black', linestyles='--')
    axe.set_aspect('equal', 'box')
    axe.set_xlabel(r'$Re(\lambda)$')
    axe.set_ylabel(r'$Im(\lambda)$')
    fig.tight_layout()
