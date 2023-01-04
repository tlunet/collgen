#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 12:07:46 2023

@author: telu
"""
import numpy as np

def matVecMul(mat, u):
    r"""
    Compute vectorized Matrix Vector Multiplication :math:`Ax` (A @ x)

    Parameters
    ----------
    mat : np.ndarray, size (nDOF, M, M) or (M, M)
        Matrix or array of matrices.
    u : np.ndarray, size (nDOF, M) or (M,)
        Vector or array of vectors.

    Returns
    -------
    out : np.ndarray, size (nDOF, M) or (M,)
        The computed matrix-vector product(s)

    Notes
    -----
    - matVecMul((nDOF, M, M), (nDOF, M)) -> (nDOF, M) <=> (M, M) @ (M,) for each nDOF
    - matVecMul((M, M), (nDOF, M)) -> (nDOF, M) <=> (M, M) @ (M,) for each nDOF
    - matVecMul((M, M), (M,)) -> (M,) <=> (M, M) @ (M,)
    """
    return np.matmul(mat, u[..., None]).squeeze(axis=-1)


def matVecInv(mat, u):
    r"""
    Compute vectorized Matrix Vector Inversion :math:`A^{-1}x` (A / x)

    Parameters
    ----------
    mat : np.ndarray, size (nDOF, M, M) or (M, M)
        Matrix or array of matrices.
    u : np.ndarray, size (nDOF, M) or (M,)
        Vector or array of vectors.

    Returns
    -------
    out : np.ndarray, size (nDOF, M) or (M,)
        The computed matrix-vector inversion(s)

    Notes
    -----
    - matVecInv((nDOF, M, M), (nDOF, M)) -> (nDOF, M) <=> (M, M) \ (M,) for each nDOF
    - matVecInv((M, M), (nDOF, M)) -> (nDOF, M) <=> (M, M) \ (M,) for each nDOF
    - matVecInv((M, M), (M,)) -> (M,) <=> (M, M) \ (M,)
    """
    try:
        return np.linalg.solve(mat, u)
    except ValueError:
        return np.linalg.solve(mat[None, ...], u)
