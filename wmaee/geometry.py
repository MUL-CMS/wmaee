
import operator
import itertools
import numpy as np
from ase import Atoms
from typing import Tuple, Dict


def find_neighbors(atoms: Atoms, rcut: float = 3.0, mindist: float = 0.1, properties: Tuple[str, ...] = ("indices", "distances", "images", "vecs")) -> Dict[int, Tuple[np.ndarray, ...]]:
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
