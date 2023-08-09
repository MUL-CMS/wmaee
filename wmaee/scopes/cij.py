"""
Convenient functions for working with elasticity.
"""

import numpy as np
from typing import *
from ase import Atoms
from typing import Tuple
from itertools import product
from frozendict import frozendict as fdict
from numpy.typing import NDArray, ArrayLike


VOIGT_INDEX_MAP: Mapping[int, Tuple[int, int]] = fdict({
        0: (0, 0),
        1: (1, 1),
        2: (2, 2),
        3: (1, 2),
        4: (0, 2),
        5: (0, 1)
})


def get_ULICS(max_eps: float = 1.5e-2) -> NDArray:
    """
    Returns ULICS (universal linear-independent cou-pling strains) used
    for deformations in the stress-strain method.
    For details see doi:10.1016/j.cpc.2009.11.017.
    :param max_eps: (float) magnitide of the stress (largest component
        in the Voigt's notation). Defaults to 1.5e-2.
    :return: np.array 6x6 matrix of ULICS
    """
    ULICS = max_eps / 6.0 * np.array([
        [1, -2, 3, -4, 5, -6],
        [2, 1, -5, -6, 4, 3],
        [3, 4, -1, 5, 6, -2],
        [4, -3 ,6, 1, -2, 5],
        [5, 6, 2, -3, -1, -4],
        [6, -5, -4, 2, -3, 1]
        ])
    return ULICS


def apply_strain(atoms: Atoms, strain: ArrayLike, div_two: bool = True) -> Atoms:
    """Applies strain to a structure

    :param atoms: initial structure
    :type atoms: ase.Atoms
    :param strain: strain to be applied
    :type strain: ArrayLike
    :param div_two: whether apply factor of 2 when converting from
        Voigt's notation to tensorial 3x3 strain. (default is True)
    :type div_two: bool
    :return: deformed structure
    :rtype: ase.Atoms
    """
    strain = np.array(strain)
    # try to convert the vector
    if strain.shape != (3, 3):
        strain = from_voigt(strain, div_two=div_two)
    deformation_matrix = np.eye(3) + strain

    new_cell = np.dot(deformation_matrix, atoms.cell.array)
    return Atoms(cell=new_cell, scaled_positions=atoms.get_scaled_positions(), numbers=atoms.numbers)


def index_from_voigt(i: int) -> Optional[Tuple[int, int]]:
    """
    Converts Voigt's (1 index) and tensorial (2 indices) notation

    :param i:  index in Voigt's notation (0..5)
    :type i: int
    :return: pair of indices (0..2)
    :rtype: Optional[Tuple[int, int]]
    """

    return VOIGT_INDEX_MAP.get(i)


def index_to_voigt(i: int, j: int) -> Optional[int]:
    """
    Converts tensorial (2 indices) to Voigt's (1 index) notation

    :param i: index in tensorial notation (0..2)
    :type i: int
    :param j: index in tensorial notation (0..2)
    :type j: int
    :return: index in Voigt's notation (0..6)
    :rtype: int
    """
    return {v: k for k, v in VOIGT_INDEX_MAP.items()}.get((i, j))


def transform_tensor(tensor: NDArray, arr: NDArray, fact_two: bool = False) -> NDArray:
    """
    Transforms tensor according to a given transformation matrix
    :param tensor: tensor to be transformed
    :type tensor: NDArray
    :param arr: transformation matrix
    :type arr: NDArray
    :param fact_two: whether to apply factors 2, 4, ... during conversion to/from Voigt's notation
    :type fact_two: bool
    :return: transformed matrix
    :rtype: NDArray
    """
    if tensor.shape in {(6, 6), (3, 3, 3, 3)}:
        if tensor.shape == (6, 6):
            voigt = True
            tensor = from_voigt(tensor, fact_two)
        else:
            voigt = False
        e = np.zeros(81).reshape((3, 3, 3 ,3))
        for i, j, k, l in product(np.arange(3), np.arange(3), np.arange(3), np.arange(3)):
            for a, b, c, d in product(np.arange(3), np.arange(3), np.arange(3), np.arange(3)):
                e[i, j, k, l] += np.around(arr[i, a]*arr[j, b]*arr[k, c]*arr[l, d]*tensor[a, b, c, d], 2)
        if voigt:
            e = to_voigt(e, fact_two)
        return e
    else:
        raise ValueError('Unknown shape of the input tensor.')


def from_voigt(m: NDArray, div_two: bool = True) -> NDArray:
    """
    Transform the strain stress array from Voigt's notation into a (3,3) matrix. In case it is a (6,6) matrix
    the function will produce a forth rank tensor.

    :param m: the stress or strain array
    :type m: NDArray
    :param div_two: divide off-diagonal elements by two (default is True)
    :type div_two: bool
    :return: a transformed tensor
    :rtype: NDArray
    """

    if m.shape == (6,):
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


def to_voigt(m: NDArray, times_two: bool = True) -> NDArray:
    """
    Transform a strain stress array tensor into Voigt's notation (3,3) -> (6,). In case {m} is a forth rank tensor
    a 2D matrix wil be produced. (3,3,3,3) -> (6,6).

    :param m: the stress or strain array
    :type m: NDArray
    :param times_two: multiply off-diagonal elements by two (default is True)
    :type times_two: bool
    :return: a transformed tensor
    :rtype: NDArray
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
           for i in range(3, 6):
               fact[i] = 2.0
        voigt = np.zeros(36).reshape((6, 6))
        for i, j in product(range(6), range(6)):
            # index_from_voigt returns an Optional type, therefore if a "None" would occur the function breaks at the unpacking step
            a, b = index_from_voigt(i)
            c, d = index_from_voigt(j)
            voigt[i, j] = m[a, b, c, d]*fact[i]*fact[j]
        return voigt
    else:
        raise ValueError('A stress matrix must be of shape (3,3)')


def project_cubic(cij: NDArray) -> NDArray:
    """
    Computes to Cij tensor projected to cubic symmetry

    :param cij: array of shape (6,6) the raw elasticity tensor
    :type cij: NDArray
    :return: the projected Cij tensor
    :rtype: NDArray
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


def project_hexagonal(cij: NDArray) -> NDArray:
    """
    Computes Cij tensor projected to hexagonal symmetry
    :param cij: of shape (6,6) the raw elasticity tensor
    :type cij: NDArray
    :return: the projected Cij tensor
    :rtype: NDArray
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
