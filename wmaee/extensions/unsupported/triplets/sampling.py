from ase import Atoms
import numpy as np
import pickle
from os.path import exists
from tqdm import tqdm
from typing import Union, Optional, Collection, Tuple, Iterable
from wmaee.core.types import Atoms

ArrayLike = Union[np.ndarray, Collection]


def sampling_function(x: ArrayLike, mu: float, num: int, pot: Tuple[float],
                      indices: Optional[bool] = False) -> ArrayLike:
    """
    Calculates the asymmetric sampling with Mie-potential weights, around the equilibrium mu and returns the samples.
    :param x: (np.ndarray) the axis whish should be sampled
    :param mu: (float) the expected minimum
    :param num: (int) the number of samples
    :param pot: (tuple of int) the powers of the Mie.Potential weighting function
    :param indices: (bool) wether the samples for the indicies of the samples should be returned
    :return: (np.ndarray or list) the samples between np.amin(x) and np.amax(x) including the boundaries
    """
    # LJ sampling
    sampling_importance = np.power(mu * np.reciprocal(x), pot[0]) - 2 * np.power(mu * np.reciprocal(x), pot[1])
    sampling_importance[sampling_importance > 0] = 0.0
    sampling_importance = np.abs(sampling_importance)
    # cute the negative values
    if np.any(np.isclose(sampling_importance, 0.0)):
        max_zero_index = np.amax(np.array(np.argwhere(np.isclose(sampling_importance, 0)).flat))
        sampling_importance = sampling_importance[max_zero_index:]
    sampling_importance /= np.amax(sampling_importance)
    # calculate the step size
    dx = (np.amax(x) - np.amin(x)) / len(x)
    integration = np.cumsum(sampling_importance * dx)
    norm = np.amax(integration)
    areas = integration / norm * num
    si = [0]
    # Find the sample indices by minimizing the distance to a certain step
    for sample in range(1, num - 1):
        re = np.abs(areas - sample)
        si.append(list(np.argwhere(re == np.amin(re)).flat)[0])
    si.append(len(areas) - 1)
    fine_axis = np.linspace(np.amin(x), np.amax(x), num=np.amax(si) + 1)
    samples = [fine_axis[fis] for fis in si]
    return samples if not indices else np.array(si)


def sample(start: float, cut: float, mu: float, num: Optional[int] = 50, pot: Tuple[float] = (5, 1)) -> np.ndarray:
    """
    Samples a half circle of triplets, from a minimum distance to the a maximum distance of "cut" and denser
    sampling around "mu", where "pot" are the weights of a Mie-Type weighting function
    :param start: (float) the minimum distance between the atoms
    :param cut: (float) the cutoff radius
    :param mu: (float) the expected minimum
    :param num: (int) the number of samples on each axis
    :param pot: (tuple of int) the powers of the Mie.Potential weighting function
    :return: (np.ndarray) a array with num**3 rows, and 3 columns which are r_ik, rjk and theta_ijk
    """
    x = np.linspace(start, cut, num=10000)
    x[np.isclose(x, 0.0)] = 1e-10
    samples = sampling_function(x, mu, num, pot)

    # we have to calculate the minimal angles needed to keep the minimum distance
    min_angles = 2 * np.arcsin(start / (2 * np.array(samples)))
    max_angles = np.zeros_like(min_angles) + np.pi
    coords = []
    num_samples = 0
    for rsample, min_angle, max_angle, in tqdm(zip(samples, min_angles, max_angles), total=num):
        angles = np.linspace(min_angle, max_angle, num=10000)
        second = np.array([rsample * np.ones_like(angles), np.zeros_like(angles), np.zeros_like(angles)])
        for ksample in samples:
            third = np.array([ksample * np.cos(angles), ksample * np.sin(angles), np.zeros_like(angles)])
            dists = np.linalg.norm(third - second, axis=0)
            if np.amax(dists) < mu:
                for ang in np.linspace(min_angle, max_angle, num=num):
                    coords.append((rsample, ksample, ang))
            else:
                si = sampling_function(dists, mu, num, pot, indices=True)
                fine_axis = np.linspace(min_angle, max_angle, num=np.amax(si) + 1)
                ang_samples = [fine_axis[fis] for fis in si]
                for ang in (ang_samples if num_samples % 2 == 0 else reversed(ang_samples)):
                    coords.append((rsample, ksample, ang))
                num_samples += 1
    return np.array(coords).T


def construct_triangle(rij: ArrayLike, rik: ArrayLike, theta: ArrayLike, origin: Tuple[float] = (0.0, 0.0, 0.0)) -> \
        Tuple[np.ndarray]:
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
              pot: Tuple[float] = (5, 1),
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
    grid = sample(start, cut, mu, num=num, pot=pot)
    triangs = np.array(construct_triangle(*grid))
    _, _, three = triangs
    cell = (np.eye(3) * (np.amax(three, axis=1) - np.amin(three, axis=1)) + np.array([vacuum] * 3)) * np.eye(3)
    pickle_file = 'grid.pickle'
    if not exists(pickle_file):
        with open(pickle_file, 'wb') as pickle_handle:
            pickle.dump(grid, pickle_handle)
    for i in range(triangs.shape[-1]):
        yield Atoms(cell=cell, symbols=species, positions=triangs[:, :, i])
