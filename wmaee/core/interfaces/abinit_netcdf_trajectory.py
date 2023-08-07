
import numpy as np
from typing import List
from ase.atoms import Atoms
from netCDF4 import Dataset
from ase.units import Hartree, Bohr
from ase.calculators.abinit import Abinit
from ase.calculators.calculator import Calculator


def make_history_filename(abinit: Abinit):
    return f"{abinit.label}o_HIST.nc"


def read_num_structures(dset: Dataset) -> int:
    n, *_ = dset.variables["ekin"].shape
    return n


def read_species(dset: Dataset) -> List[int]:
    zvals = [int(ordinal) for ordinal in dset.variables["znucl"][:]]
    return [zvals[int(typid) - 1] for typid in dset.variables["typat"][:]]


def read_cell(dset: Dataset, index: int) -> np.ndarray:
    return np.array(dset.variables["rprimd"][index, :, :]) * Bohr


def read_positions(dset: Dataset, index: int) -> np.ndarray:
    return np.array(dset.variables["xcart"][index, :, :]) * Bohr


class AbinitNetCDFDummyCalculator(Calculator):

    implemented_properties = ['energy', 'forces', 'stress', 'kinetic_energy']

    def __init__(self, dset: Dataset, index: int, **kwargs):
        super().__init__(**kwargs)
        self._dset = dset
        self._index = index

    def get_forces(self, atoms=None):
        return np.array(self._dset.variables["fcart"][self._index, :, :]) * Hartree / Bohr

    def get_stress(self, atoms=None):
        return np.array(self._dset.variables["strten"][self._index, :]) * Hartree / (Bohr**3)

    def get_potential_energy(self, atoms=None, force_consistent=False):
        return self.get_total_energy(atoms=atoms) - self.get_kinetic_energy(atoms=atoms)

    def get_kinetic_energy(self, atoms=None):
        return self._dset.variables["ekin"][self._index] * Hartree

    def get_total_energy(self, atoms=None):
        return self._dset.variables["etotal"][self._index] * Hartree


def read_abinit_netcdf_trajectory(abinit: Abinit):
    dset = Dataset(make_history_filename(abinit))
    numbers = read_species(dset)
    for i in range(read_num_structures(dset)):
        atoms = Atoms(cell=read_cell(dset, i), positions=read_positions(dset, i), numbers=numbers, calculator=AbinitNetCDFDummyCalculator(dset, i))
        yield atoms
    dset.close()