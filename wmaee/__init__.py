# This is to hide moudle data
# plotting functions
from .wmaee import plot_projected_dos, plot_total_dos
from plotly.graph_objects import Figure, Scatter, Contour
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

# Forward system functions, since it should be for dummies
from os import getcwd as current_directory, listdir, walk, mkdir
from os.path import isfile, isdir, join, exists

# Greet from
from .wmaee import greet as _greet

_greet()