
import numpy as np
from ase.atoms import Atoms
from netCDF4 import Dataset
from typing import List, Optional
from ase.units import Hartree, Bohr
from ase.calculators.abinit import Abinit
from ase.calculators.calculator import Calculator


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


class AbinitNetCDF4TrajectoryReader:

    def __init__(self, abinit: Abinit):
        self._abinit: Abinit = abinit
        self._dset: Optional[Dataset] = None
        self._numbers: Optional[List[int]] = None
        self._num_structures: Optional[int] = None
        self._indexer: Optional[List[int]] = None
        self._open(f"{self._abinit.label}o_HIST.nc")

    def _open(self, filename: str):
        self._dset = Dataset(filename)
        self._numbers = self.read_species()
        self._num_structures = self.read_num_structures()
        self._indexer = list(range(self._num_structures))

    def read_species(self) -> List[int]:
        zvals = [int(ordinal) for ordinal in self._dset.variables["znucl"][:]]
        return [zvals[int(typid) - 1] for typid in self._dset.variables["typat"][:]]

    def read_num_structures(self) -> int:
        n, *_ = self._dset.variables["ekin"].shape
        return n

    def read_cell(self, index: int) -> np.ndarray:
        return np.array(self._dset.variables["rprimd"][index, :, :]) * Bohr

    def read_positions(self, index: int) -> np.ndarray:
        return np.array(self._dset.variables["xcart"][index, :, :]) * Bohr

    def read_atoms(self, index: int) -> Atoms:
        return Atoms(
            cell=self.read_cell(index),
            positions=self.read_positions(index),
            numbers=self._numbers.copy(),
            calculator=AbinitNetCDFDummyCalculator(self._dset, index)
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._dset.isopen():
            self._dset.close()

    def __len__(self):
        return self._num_structures

    def __getitem__(self, item):
        if isinstance(item, slice):
            return tuple(map(self.read_atoms, self._indexer[item]))
        else:
            return self.read_atoms(item)

    def __iter__(self):
        for i in range(self._num_structures):
            yield self[i]
