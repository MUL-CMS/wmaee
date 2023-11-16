import operator
import itertools
import numpy as np
from ase import Atoms
from typing import Union, Tuple, Dict, Any

from wmaee.core.config import is_pmg_avail
if is_pmg_avail():
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor

def find_neighbors(atoms: Union[Atoms, Any],
                   rcut: float = 3.0,
                   mindist: float = 0.1,
                   properties: Tuple[str, ...] = ('indices', 'distances', 'images', 'vecs')
                   ) -> Dict[int, Tuple[np.ndarray, ...]]:
    """
    Find neighbors for each atom within a specified cutoff distance and minimum distance.
    
    This function applies the minimum image convention to consider periodic boundary conditions.

    Parameters
    ----------
    atoms : Union[Atoms, Any]
        The atoms object to query. It can be an Atoms object or any other compatible type.
    rcut : float, optional
        The cutoff distance for neighbor search, by default 3.0.
    mindist : float, optional
        The minimum distance for neighbor search, by default 0.1.
    properties : Tuple[str, ...], optional
        An iterable of properties to return for each atom. 
        Available options are 'indices', 'distances', 'images', 'vecs', by default ('indices', 'distances', 'images', 'vecs').

    Returns
    -------
    Dict[int, Tuple[np.ndarray, ...]]
        A dictionary where the keys represent the atom ids, and the values are tuples containing the desired properties.
        - 'indices': Indices of neighboring atoms.
        - 'distances': Distances to neighboring atoms.
        - 'images': Indices of periodic images of neighboring atoms.
        - 'vecs': Vectors pointing from the central atom to its neighbors.
    """
    if is_pmg_avail() and isinstance(atoms, Structure):
        atoms = AseAtomsAdaptor.get_atoms(atoms)

    a, b, c = atoms.cell
    stride = len(atoms)
    image_vector = np.zeros((stride * 27, 3), dtype=float)
    indices = np.hstack([np.arange(len(atoms))] * 27)
    image_indices = np.zeros_like(image_vector, dtype=int)

    for i, (ta, tb, tc) in enumerate(itertools.product(*[[-1, 0, 1]] * 3)):
        start, end = (i * stride), ((i + 1) * stride)
        image_vector[start:end] = a * ta + b * tb + c * tc
        image_vector[start:end] += atoms.positions
        image_indices[start:end] = [ta, tb, tc]

    result = dict()
    position_vector = np.zeros_like(image_vector)

    for i, pos in enumerate(atoms.positions):
        position_vector[:] = pos
        distance_vectors = position_vector - image_vector
        distances = np.linalg.norm(distance_vectors, axis=1)
        mask = np.logical_and(distances >= mindist, distances < rcut)
        result[i] = dict(indices=indices[mask], distances=distances[mask], images=image_indices[mask], vecs=distance_vectors[mask, :])

    mapper = operator.itemgetter(*properties)
    return {i: mapper(d) for i, d in result.items()}