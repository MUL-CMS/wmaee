import numpy as np
from wmaee.extensions.unsupported.triplets.sampling import sampling_function, ArrayLike
from wmaee.core.types import Atoms, Collection, Union, Optional, Tuple, Iterable


def sample_pair_axes(start: float, cut: float, mu: float, num: Optional[int] = 50,
                     pot: Tuple[float, float] = (5, 1)) -> ArrayLike:
    """
    Samples a half circle of triplets, from a minimum distance to the a maximum distance of "cut" and denser
    sampling around "mu", where "pot" are the weights of a Mie-Type weighting function
    :param start: (float) the minimum distance between the atoms
    :param cut: (float) the cutoff radius
    :param mu: (float) the expected minimum
    :param num: (int) the number of samples on each axis
    :param pot: (tuple of int) the powers of the Mie.Potential weighting function
    :return: (np.ndarray) sample of the axis
    """
    x = np.linspace(start, cut, num=10000)
    x[np.isclose(x, 0.0)] = 1e-10
    samples = sampling_function(x, mu, num, pot)
    return samples


def pairs(start: float, cut: float, mu: float, species: Union[Collection[str], Tuple[str]], num: Optional[int] = 50,
          pot: Tuple[float, float] = (5, 1), vacuum: Optional[float] = 20) -> Iterable[Atoms]:
    samples = sample_pair_axes(start, cut, mu, num=num, pot=pot)

    coords = None
    return samples
