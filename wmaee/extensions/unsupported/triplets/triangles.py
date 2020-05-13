import pickle
import numpy as np
from ase import Atoms as AseAtoms
from os.path import exists
from wmaee.extensions.unsupported.triplets.sampling import sample_triangle, ArrayLike
from wmaee.core.types import Union, Collection, Optional, Iterable, Tuple, Atoms


def construct_triangle(rij: ArrayLike, rik: ArrayLike, theta: ArrayLike,
                       origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)) -> Tuple[np.ndarray]:
    """
    Constructs triangles, where first atom sits at origin, the second one along the x-axis and the third one
    is sampled around in a half circle
    :param rij: (np.ndarray) distance between atom 1 and 2
    :param rik: (np.ndarray) distance between atom 1 and 3
    :param theta: (np.ndarray) angle which is opposite of the hypothenuse
    :param origin: (tuple of float) a shift of the origin
    :return: (tuple np.ndarray) cordinates of atom 1, 2, and three
    """
    origin = np.array(origin)
    first = np.zeros((3, len(rij)))
    second = np.zeros((3, len(rij)))
    third = np.zeros((3, len(rij)))
    second[0, :] = rij
    third[0, :] = rik * np.cos(theta)
    third[1, :] = rik * np.sin(theta)
    for i, coord in enumerate(origin):
        first[i, :] += coord
        second[i, :] += coord
        third[i, :] += coord
    return first, second, third


def triangles(start: float, cut: float, mu: float, species: Union[Collection[str], Tuple[str]], num: Optional[int] = 50,
              pot: Tuple[float, float] = (5, 1),
              vacuum: Optional[float] = 20) -> Iterable[Atoms]:
    """
    :param species: (list of str) the species strings
    :param start: (float) the minimum distance between the atoms
    :param cut: (float) the cutoff radius
    :param mu: (float) the expected minimum
    :param num: (int) the number of samples on each axis
    :param pot: (tuple of int) the powers of the Mie.Potential weighting function
    :param vacuum: (float) how much vacuum will be added to pad from the next periodic image
    :return:
    """
    grid = sample_triangle(start, cut, mu, num=num, pot=pot)
    triangs = np.array(construct_triangle(*grid))
    _, _, three = triangs
    cell = (np.eye(3) * (np.amax(three, axis=1) - np.amin(three, axis=1)) + np.array([vacuum] * 3)) * np.eye(3)
    pickle_file = 'grid.pickle'
    if not exists(pickle_file):
        with open(pickle_file, 'wb') as pickle_handle:
            pickle.dump(grid, pickle_handle)
    for i in range(triangs.shape[-1]):
        yield AseAtoms(cell=cell, symbols=species, positions=triangs[:, :, i])
