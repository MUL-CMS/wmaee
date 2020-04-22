
from pymatgen.io.ase import AseAtomsAdaptor
from collections.abc import Iterable

def pymatgen_to_ase(structure):
    """
    Converts a `pymatgen.Structure` object to a `ase.Atoms` object
    :param structure: (pymatgen.Structure) the structure to convert
    :return: (ase.Atom) the converted object
    """
    return AseAtomsAdaptor.get_atoms(structure)

def ase_to_pymatgen(atoms):
    """
     Converts  a `ase.Atoms` object to a `pymatgen.Structure` object
    :param atoms: (ase.Atom) the structure to convert
    :return: (pymatgen.Structure) the converted object
    """
    return AseAtomsAdaptor.get_structure(atoms)

def collection(a):
    return [a] if isinstance(a, str) else (a if isinstance(a, (set, list, tuple)) else [a])

def unpack_single(a):
    return a[0] if len(a) == 1 else a