from pyiron.atomistics.structure.atoms import Atoms as IronAtoms
from ase import Atoms as AseAtoms
from wmaee.core.common import working_directory
from pymatgen import Structure as PymatgenStructure
from typing import Union, Dict, Tuple, Collection, NoReturn, Optional, List, Iterator, TextIO, Callable

Atoms = Union[PymatgenStructure, IronAtoms, AseAtoms]
Directory = Union[None, str, working_directory]