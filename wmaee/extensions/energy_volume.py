
"""
Convenient functions for energy-volume calculations for the WMaE exercises.
"""

import os
import numpy as np

from ase import Atoms
from os.path import join, exists
from functools import partial
from numpy.typing import ArrayLike
from scipy.optimize import curve_fit


def energy_volume_curve(name: str, atoms: Atoms, deltas=None, incar=None, kpts=(1, 1, 1), show_output=False, ncpus=2,
                              xc="PBE"):
    pass
    # """
    # Performs a E-V curve calculation by loading a cell with hydrostatic strain. Returns a list of volumes and energies
    # in A^3 and eV
    #
    # :param name: (str) the name if the folder where to locate the data
    # :param atoms: (pymatgen.Structure or ase.Atoms) the structure
    # :param deltas: (list, tuple or np.array) the magnitude of the hydrostatic strains (default: linspace(-0.05, 0.05,10))
    # :param delta: (pymatgen.io.vasp.Kpoints) a KPOINTS definition to use (default: Kpoints.automatic_density_by_vol(32000))
    # :param show_output: (bool) wether to show VASP output or not
    # :param ncpus: (int) number of cpus to use (default is 2)
    # :param xc: (str) the XC type identifier (default is "pbe")
    # :return: (list) the calculated data
    # """
    #
    # finished = lambda : exists('vasprun.xml')
    # if incar is None:
    #     incar = Incar(dict(
    #         NSW=0,
    #         LREAL='Auto',
    #         LCHARG='False',
    #         LWAVE='False'
    #     ))
    #
    # if kpoints is None:
    #     kpoints = Kpoints.automatic_density_by_vol(32000)
    #
    # if deltas is None:
    #     deltas = np.linspace(-0.05, 0.05, 10)
    #
    # hydrostatic_strain = lambda delta: np.eye(3) * delta
    # folder_name = lambda seq: 'strain-%i' % seq
    #
    # with folder(name):
    #     data = []
    #     for i, delta in enumerate(tqdm(deltas, desc='Progress')):
    #         with folder(folder_name(i)):
    #             # speed up the calculation by symliniking WAVECAR
    #             if i == 0:
    #                 incar['LWAVE'] = True
    #             else:
    #                 incar['LWAVE'] = False
    #                 # make symlink
    #                 os.system('ln -s %s %s' % (join('..', folder_name(0), 'WAVECAR'), 'WAVECAR'))
    #             deformed = apply_strain(structure, hydrostatic_strain(delta))
    #             # print(deformed.lattice.matrix/4.05)
    #             energy = full_run(VASPInput(incar, deformed, kpoints, xc_func=xc_func), cpus=cpus,
    #                               show_output=show_output).final_energy if not finished() else parse_output().final_energy
    #             data.append((deformed.lattice.volume, energy))
    #     return list(zip(*data))


def birch_murnaghan(volumes: ArrayLike, e0: float, v0: float, b0: float, bp: float):
    two_thirds = 2/3
    return e0 + (9 * v0 * b0 / 16) * (bp * ((volumes / v0) ** two_thirds - 1) ** 3 + ((volumes / v0) ** two_thirds - 1) ** 2 * (6 - 4 * (volumes / v0) ** two_thirds))


def birch_murnaghan_fit(volumes, energies, p0=None):
    """
    Performs a Birch-Murnaghan-fit on a given set of volumes and energies
    :param volumes: volumes
    :type volumes: ArrayLike
    :param energies: energies
    :type energies: ArrayLike
    :param p0: (iterable) starting parameters in the order E0_guess, V0_guess, B0_guess, Bp_guess (default: None)
    :param return_analytic: (bool) wether also to return a SymPy analytic expression of the fit
    :returns: (optimum:np.array, E(V):callable, optional analytic:Sympy.Expression)
    """

    if p0 is None:
        V0_guess = np.amin(volumes)
        E0_guess = np.amin(energies)
        B0_guess = 100
        Bp_guess = 1
        p0 = (E0_guess, V0_guess, B0_guess, Bp_guess)

    optimum, residuals = curve_fit(birch_murnaghan, np.array(volumes), np.array(energies), p0=p0)
    return optimum, lambda volumes: birch_murnaghan(volumes, *optimum)