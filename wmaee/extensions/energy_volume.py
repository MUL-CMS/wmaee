import os
import numpy as np
from wmaee.core.common import tqdm, working_directory as folder
from wmaee.vasp import VASPInput, parse_output, full_run
from wmaee.utils import ase_to_pymatgen
from wmaee.extensions.elasticity import apply_strain
from pymatgen.io.vasp import Incar, Kpoints
from ase import Atoms
from os.path import join, exists
from functools import partial
from scipy.optimize import curve_fit
from sympy import symbols, lambdify, S

def energy_volume_curve(name, structure, deltas=None, incar=None, kpoints=None, show_output=False, cpus=2,
                             xc_func='gga'):
    """
    Performs a E-V curve calculation by loading a cell with hydrostatic strain. Returns a list of volumes and energies
    in A^3 and eV

    :param name: (str) the name if the folder where to locate the data
    :param structure: (pymatgen.Structure or ase.Atoms) the structure
    :param deltas: (list, tuple or np.array) the magnitude of the hydrostatic strains (default: linspace(-0.05, 0.05,10))
    :param delta: (pymatgen.io.vasp.Kpoints) a KPOINTS definition to use (default: Kpoints.automatic_density_by_vol(32000))
    :param incar: (pymatgen.io.vasp.Incar) a INCAR definition to use
    :param show_output: (bool) wether to show VASP output or not
    :param cpus: (int) number of cpus to use (default: 2)
    :param xc_func: (str) the XC type identifier (default: gga)
    :return: (list) the calculated data
    """
    if isinstance(structure, Atoms):
        structure = ase_to_pymatgen(structure)

    finished = lambda : exists('vasprun.xml')
    if incar is None:
        incar = Incar(dict(
            NSW=0,
            LREAL='Auto',
            LCHARG='False',
            LWAVE='False'
        ))

    if kpoints is None:
        kpoints = Kpoints.automatic_density_by_vol(32000)

    if deltas is None:
        deltas = np.linspace(-0.05, 0.05, 10)

    hydrostatic_strain = lambda delta: np.eye(3) * delta
    folder_name = lambda seq: 'strain-%i' % seq

    with folder(name):
        data = []
        for i, delta in enumerate(tqdm(deltas, desc='Progress')):
            with folder(folder_name(i)):
                # speed up the calculation by symliniking WAVECAR
                if i == 0:
                    incar['LWAVE'] = True
                else:
                    incar['LWAVE'] = False
                    # make symlink
                    os.system('ln -s %s %s' % (join('..', folder_name(0), 'WAVECAR'), 'WAVECAR'))
                deformed = apply_strain(structure, hydrostatic_strain(delta))
                # print(deformed.lattice.matrix/4.05)
                energy = full_run(VASPInput(incar, deformed, kpoints, xc_func=xc_func), cpus=cpus,
                                  show_output=show_output).final_energy if not finished() else parse_output().final_energy
                data.append((deformed.lattice.volume, energy))
        return list(zip(*data))

def birch_murnaghan_fit(volumes, energies, p0=None, return_analytic=False):
    """
    Performs a Birch-Murnaghan-fit on a given set of volues and energies
    :param volumes: (iterable) volumes
    :param energies (iterable) energies
    :param p0: (iterable) starting parameters in the order E0_guess, V0_guess, B0_guess, Bp_guess (default: None)
    :param return_analytic: (bool) wether also to return a SymPy analytic expression of the fit
    :returns: (optimum:np.array, E(V):callable, optional analytic:Sympy.Expression)
    """

    # We define the variables at first
    E0, V0, B0, V, Bp = symbols('E_0 V_0 B_0 V B_p')
    E = E0 + (9 * V0 * B0 / 16) * (
                Bp * ((V / V0) ** S('2/3') - 1) ** 3 + ((V / V0) ** S('2/3') - 1) ** 2 * (6 - 4 * (V / V0) ** S('2/3')))

    equation_of_state = lambdify(E.free_symbols, E, dummify=False)
    if p0 is None:
        V0_guess = volumes[energies.index(min(energies))]
        E0_guess = min(energies)
        B0_guess = 100
        Bp_guess = 1
        p0 = (E0_guess, V0_guess, B0_guess, Bp_guess)
    else:
        pass

    def eos(x, E0, V0, B0, Bp):
        kwargs = dict(V=x, E_0=E0, V_0=V0, B_0=B0, B_p=Bp)
        return equation_of_state(**kwargs)

    optimum, residuals = curve_fit(eos, volumes, energies, p0=p0)
    # Unpack the results
    E0_value, V0_value, B0_value, Bp_value = optimum
    # make a function E_(V)
    E_ = partial(eos, E0=E0_value, V0=V0_value, B0=B0_value, Bp=Bp_value)
    # before fitting we have to do a small trick
    result = [optimum, E_]
    if return_analytic:
        result.append(E.subs(dict(E_0=E0_value, V_0=V0_value, B_0=B0_value, B_p=Bp_value)))

    return tuple(result)