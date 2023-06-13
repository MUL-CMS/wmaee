"""
Convenient functions for working with elasticity.
"""

import numpy as np
from typing import Tuple
from itertools import product


def index_from_voigt(i: int) -> Tuple[int, int]:
    """
    Converts Voigt's (1 index) and tensorial (2 indices) notation

    Parameters
    ----------
    i : int
        Index in Voigt's notation (0..5).

    Raises
    ------
    ValueError
        Wrong value of index.

    Returns
    -------
    (int, int)
        Index of the (3x3) tensor component.

    """
    if i == 0:
        return 0, 0
    elif i == 1:
        return 1, 1
    elif i == 2:
        return 2, 2
    elif i == 3:
        return 1, 2
    elif i == 4:
        return 0, 2
    elif i == 5:
        return 0, 1
    else:
        raise ValueError('Unknown value of index')



def index_to_voigt(i: int, j: int) -> int:
    """
    Converts tensorial (2 indices) to Voigt's (1 index) notation.

    Parameters
    ----------
    i : int
        Index in tensorial notation (0..2).
    j : int
        index in tensorial notation (0..2).

    Raises
    ------
    ValueError
        Input values are not (0..2) integers.

    Returns
    -------
    int
        Index in Voigt's notation (0..6).

    """
    if (i, j) in product(np.arange(3), np.arange(3)):
        if i == j:
            return i
        elif i+j == 3:
            return 4
        elif i+j == 2:
            return 5
        elif i+j == 3:
            return 6
    else:
        raise ValueError('Unknown value of indices')


def transform_tensor(t: np.ndarray, a: np.ndarray, fact_two: bool = False) -> np.ndarray:
    """
    Transforms tensor according to a given transformation matrix

    Parameters
    ----------
    t : numpy.ndarray
        Tensor to be transformed represented eiter by (6x6) or 
        (3x3x3x3) matrix.
    a : numpy.ndarray
        Transformation matrix.
    fact_two : bool, optional
         Whether to apply factors 2, 4, ... during conversion 
         to/from Voigt's notation. The default is False.

    Raises
    ------
    ValueError
        Wrong shape of the input tensor t.

    Returns
    -------
    e : numpy.ndarray
        Transformed matrix in the same shape as the input
        tensor t.

    """
    if t.shape in [(6, 6), (3, 3, 3, 3)]:
        if t.shape == (6, 6):
            voigt = True
            t = from_voigt(t, fact_two)
        e = np.zeros(81).reshape((3, 3, 3 ,3))
        for i, j, k, l in product(np.arange(3), np.arange(3), np.arange(3), np.arange(3)):
            for a, b, c, d in product(np.arange(3), np.arange(3), np.arange(3), np.arange(3)):
                e[i, j, k, l] += np.around(a[i, a] * a[j, b] * a[k, c] * a[l, d] * t[a, b, c, d], 2)
        if voigt:
            e = to_voigt(e, fact_two)
        return e
    else:
        raise ValueError('Unknown shape of the input tensor.')


def from_voigt(m: np.ndarray, div_two: bool = True) -> np.ndarray:
    """
    TODO: add docstring here (but consitently)
    """
    if m.shape == (6, ):
        # 2nd rank tensor: vector 6x1 -> matrix 3x3
        e = m.copy()
        if div_two:
            e[3] /= 2.0
            e[4] /= 2.0
            e[5] /= 2.0
        m = np.array([
            [e[0], e[5], e[4]],
            [e[5], e[1], e[3]],
            [e[4], e[3], e[2]]
        ])
        return m
    elif m.shape == (6, 6):
        # 4nd rank tensor: matrix 6x6 -> matrix 3x3x3x3
        fact = np.ones(6)
        if div_two:
           for i in np.arange(3, 6):
               fact[i] = 0.5
        e = np.zeros(81).reshape((3, 3, 3 ,3))
        for i, j in product(np.arange(6), np.arange(6)):
            a, b = index_from_voigt(i)
            c, d = index_from_voigt(j)
            e[a, b, c, d] = m[i, j]*fact[i]*fact[j]
            e[a, b, d, c] = m[i, j]*fact[i]*fact[j]
            e[b, a, c, d] = m[i, j]*fact[i]*fact[j]
            e[b, a, d, c] = m[i, j]*fact[i]*fact[j]        
        return e
    else:
        raise ValueError('Unknown shape of the input data')    


def to_voigt(m: np.ndarray, times_two: bool = True) -> np.ndarray:
    """
    TODO: add docstring here (but consitently)
    """
    if m.shape == (3, 3):
        voigt = np.array([
            m[0, 0],
            m[1, 1],
            m[2, 2],
            (m[2, 1] + m[1, 2]) / 2.0,
            (m[2, 0] + m[0, 2]) / 2.0,
            (m[1, 0] + m[0, 1]) / 2.0
        ])
        if times_two:
            voigt[3] *= 2.0
            voigt[4] *= 2.0
            voigt[5] *= 2.0
        return voigt
    if m.shape == (3, 3, 3, 3):
        fact = np.ones(6)
        if times_two:
           for i in np.arange(3, 6):
               fact[i] = 2.0
        voigt = np.zeros(36).reshape((6, 6))
        for i, j in product(np.arange(6), np.arange(6)):
            a, b = index_from_voigt(i)
            c, d = index_from_voigt(j)
            voigt[i, j] = m[a, b, c, d]*fact[i]*fact[j] 
        return voigt
    else:
        raise ValueError('A stress matrix must be of shape (3,3)')


def project_cubic(cij: np.ndarray) -> np.ndarray:
    """
    Computes to Cij tensor projected to cubic symmetry
    :param cij: array of shape (6,6) the raw elasticity tensor
    :type cij: np.ndarrays
    :return: the projected Cij tensor
    :rtype: np.ndarray
    """
    cij = np.array(cij)
    projected_cij = np.zeros((6, 6))
    projected_cij[0, 0] = np.around(np.mean([cij[i, i] for i in np.arange(3)]), 2)
    for i in np.arange(1, 3):
        projected_cij[i, i] = projected_cij[0, 0]
    projected_cij[3, 3] = np.around(np.mean([cij[i, i] for i in np.arange(3, 6)]), 2)
    for i in np.arange(4, 6):
        projected_cij[i, i] = projected_cij[3, 3]
    projected_cij[1, 0] = np.around(np.mean([cij[1, 0], cij[2, 0], cij[2, 1]]), 2)
    projected_cij[2, 0] = projected_cij[1, 0]
    projected_cij[2, 1] = projected_cij[1, 0]
    projected_cij[0, 1] = projected_cij[1, 0]
    projected_cij[0, 2] = projected_cij[1, 0]
    projected_cij[1, 2] = projected_cij[1, 0]
    return projected_cij


def project_hexagonal(cij: np.ndarray) -> np.ndarray:
    """
    Computes Cij tensor projected to hexagonal symmetry
    :param cij: of shape (6,6) the raw elasticity tensor
    :type cij: np.ndarray
    :return: the projected Cij tensor
    :rtype: np.ndarray
    """

    # convert from Cij to cij_hat, Moakher & Norris, Eq. 7
    cij_hat = np.array(cij).copy()
    for j in np.arange(3, 6):
        for i in np.arange(3):
            cij_hat[i, j] *= np.sqrt(2)
            cij_hat[j, i] *= np.sqrt(2)
        for i in np.arange(3, 6):
            cij_hat[i, j] *= 2
    # projecting on hexagonal symmetry, Moakher & Norris, Eq. A11a, A11b
    c11st = (3*cij_hat[0,0]+3*cij_hat[1,1]+2*cij_hat[0,1]+2*cij_hat[5,5])/8.
    c66st = (cij_hat[0,0]+cij_hat[1,1]-2*cij_hat[0,1]+2*cij_hat[5,5])/4.
    # Moakher & Norris, Eq. A15         
    projected_cij = np.zeros((6, 6))
    projected_cij[0, 0] = np.around(c11st, 2)
    projected_cij[1, 1] = projected_cij[0, 0]
    projected_cij[2, 2] = np.around(cij_hat[2, 2], 2)
    projected_cij[3, 3] = np.around(np.mean([cij_hat[3, 3], cij_hat[4, 4]]), 2)
    projected_cij[4, 4] = projected_cij[3, 3]
    projected_cij[5, 5] = np.around(c66st, 2)
    projected_cij[1, 0] = np.around(c11st-c66st, 2)
    projected_cij[0, 1] = projected_cij[1, 0]
    projected_cij[2, 0] = np.around(np.mean([cij_hat[2, 0], cij_hat[2, 1]]), 2)
    projected_cij[2, 1] = projected_cij[2, 0]
    projected_cij[0, 2] = projected_cij[2, 0]
    projected_cij[1, 2] = projected_cij[2, 0]
    # convert back,Moakher & Norris, Eq. 7
    for j in range(3, 6):
        for i in range(3):
            projected_cij[i, j] /= np.sqrt(2)
            projected_cij[j, i] /= np.sqrt(2)
        for i in range(3, 6):
            projected_cij[i, j] /= 2

    return projected_cij
