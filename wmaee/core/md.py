"""
This module's purpose is to gather the functionality from ase and forward it
"""

import functools
import numpy as np
from ase import units
from ase.md.npt import NPT
from ase.io import Trajectory
from ase.md.langevin import Langevin
from ase.md.andersen import Andersen
from ase.md.verlet import VelocityVerlet
from ase.md.nptberendsen import NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from wmaee.core.interfaces.calculation import MDCalculation
from ase.optimize import QuasiNewton, LBFGS, LBFGSLineSearch
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, ZeroRotation, Stationary


def get_property(calc: MDCalculation, getter=None) -> np.ndarray:
    return np.array([getter(atoms) for atoms in calc.get_trajectory()])


def steps(calc: MDCalculation) -> np.ndarray:
    return np.array([i for i, _ in enumerate(calc.get_trajectory())])


energies = functools.partial(get_property, getter=lambda atoms: atoms.get_total_energy())

temperatures = functools.partial(get_property, getter=lambda atoms: atoms.get_temperature())

volumes = functools.partial(get_property, getter=lambda atoms: atoms.get_cell().volume)


__all__ = ["units", "Trajectory", "Langevin", "Andersen", "NVTBerendsen", "NPTBerendsen",
           "MaxwellBoltzmannDistribution", "VelocityVerlet", "ZeroRotation", "Stationary", "NPT", "QuasiNewton",
           "LBFGS", "LBFGSLineSearch", "get_property", "steps", "energies", "temperatures", "volumes"]



