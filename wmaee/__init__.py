# import nglview to initilize JavaScript FrontendComms and KernelComms
import nglview as _nglview
__ = _nglview.demo()
# plotting functions
# core and convenience functionality
from wmaee.core.common import working_directory, tqdm
from wmaee.utils import ase_to_pymatgen, pymatgen_to_ase, pyiron_to_pymatgen, pymatgen_to_pyiron, ase_to_pyiron, pyiron_to_ase, to_pyiron, to_ase, to_pymatgen
#utility functions
from wmaee.vis import view
# VASP related stuff
from wmaee.vasp import VASPInput, VASPOutput, parse_output, write_input, full_run, vasp, vasp_interactive

# LAMMPS related stuff based on pyiron
from wmaee.lammps import LAMMPSCalculation
# import the shell
from wmaee.core.runner import Shell

#forward ASE classes
from ase import Atoms
from ase.io import read, write
from ase.build import bulk

# Forward pymatgen classes
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar, Outcar, Oszicar, Vasprun
from pymatgen import Structure, Lattice, Orbital, Spin

# Forward system functions, since it should be for dummies
from os import getcwd as current_directory, listdir, walk, mkdir
from os.path import isfile, isdir, join, exists
# Greet the user
from wmaee.wmaee import greet as _greet

# register view wrapper

_greet()