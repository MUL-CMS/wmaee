import numpy as np
from ase import Atoms
from pymatgen import Lattice, Structure
from sympy import Matrix

def apply_strain(structure, strain, div_two=True):
    strain = np.array(strain)
    # try to convert the vector
    if strain.shape != (3, 3):
        strain = from_voigt(strain, div_two=div_two)
    deformation_matrix = np.eye(3) + strain

    if isinstance(structure, Atoms):
        result = structure.copy()
        result.cell = np.dot(deformation_matrix, result.cell)
        result.set_scaled_positions(structure.get_scaled_positions())
    elif isinstance(structure, Structure):
        result = structure.copy()
        matrix = result.lattice.matrix
        new_lattice = Lattice(np.dot(deformation_matrix, matrix))
        result.lattice = new_lattice
    else:
        raise TypeError()
    return result

def from_voigt(m, div_two=True):
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

def to_voigt(m, times_two=True):
    if m.shape != (3, 3):
        raise ValueError('A stress matrix must be of shape (3,3)')
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

def project_cubic(cij, pretty=True):
    """
    Computes to Cij tensor projected to cubic symmetry
    :param cij: (arraylike) of shape (6,6) the raw elasticity tensor
    :param pretty: (bool) if the Cij tensor should be wrapped into a sympy.Matrix for nicer view in Jupyter notebooks
    :return: (numpy.ndarray or sympy.Matrix) the projected Cij tensor
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
    return Matrix(projected_cij) if pretty else projected_cij

def project_hexagonal(cij, pretty=True):
    """
    Computes Cij tensor projected to hexagonal symmetry
    :param cij: (arraylike) of shape (6,6) the raw elasticity tensor
    :param pretty: (bool) if the Cij tensor should be wrapped into a sympy.Matrix for nicer view in Jupyter notebooks
    :return: (numpy.ndarray or sympy.Matrix) the projected Cij tensor
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
    return Matrix(projected_cij) if pretty else projected_cij