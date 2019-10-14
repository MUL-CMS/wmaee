# This is to hide moudle data
# plotting functions
from .wmaee import plot_projected_dos, plot_total_dos
#utility functions
from .wmaee import view, ase_to_pymatgen, pymatgen_to_ase, working_directory

# VASP related stuff
from .wmaee import full_run, vasp, write_input, VASPInput, parse_output

#forward ASE classes
from ase import Atoms
from ase.io import read, write
from ase.build import bulk

# Forward pymatgen classes
from pymatgen import Structure, Lattice, Orbital, Spin
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar, Outcar, Oszicar, Vasprun


# Greet from
from .wmaee import greet as _greet

_greet()