from pymatgen.io.ase import AseAtomsAdaptor
from pyiron.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, ovito_to_pyiron, pyiron_to_ovito, \
    pyiron_to_pymatgen, pymatgen_to_pyiron, Atoms as IronAtoms
from typing import Any, Collection
from pymatgen import Structure
from ase import Atoms as AseAtoms
from typing import Union, Callable
from functools import wraps


def to_ase(s: Union[Structure, AseAtoms, IronAtoms]) -> AseAtoms:
    """
    Takes a possible structure type and returns an ase.Atoms object
    :param s: (Structure, ase.Atoms or pyiron.Atoms) the structure type
    :return: the ase.Atoms object
    """
    if isinstance(s, IronAtoms):
        return pyiron_to_ase(s)
    elif isinstance(s, Structure):
        return pymatgen_to_ase(s)
    elif isinstance(s, AseAtoms):
        return s
    else:
        raise TypeError

def to_pymatgen(s: Union[Structure, AseAtoms, IronAtoms]) -> Structure:
    """
    Takes a possible structure type and returns an Structure object
    :param s: (Structure, ase.Atoms or pyiron.Atoms) the structure type
    :return: the Structure object
    """
    if isinstance(s, IronAtoms):
        return pyiron_to_pymatgen(s)
    elif isinstance(s, Structure):
        return s
    elif isinstance(s, AseAtoms):
        return ase_to_pymatgen(s)
    else:
        raise TypeError

def to_pyiron(s: Union[Structure, AseAtoms, IronAtoms]) -> IronAtoms:
    """
    Takes a possible structure type and returns an pyiron.Atoms object
    :param s: (Structure, ase.Atoms or pyiron.Atoms) the structure type
    :return: the pyiron.Atoms object
    """
    if isinstance(s, IronAtoms):
        return s
    elif isinstance(s, Structure):
        return pymatgen_to_pyiron(s)
    elif isinstance(s, AseAtoms):
        return ase_to_pyiron(s)
    else:
        raise TypeError

def pymatgen_to_ase(structure: Structure) -> AseAtoms:
    """
    Converts a `pymatgen.Structure` object to a `ase.Atoms` object
    :param structure: (pymatgen.Structure) the structure to convert
    :return: (ase.Atom) the converted object
    """
    return AseAtomsAdaptor.get_atoms(structure)


def ase_to_pymatgen(atoms: AseAtoms) -> Structure:
    """
     Converts  a `ase.Atoms` object to a `pymatgen.Structure` object
    :param atoms: (ase.Atom) the structure to convert
    :return: (pymatgen.Structure) the converted object
    """
    return AseAtomsAdaptor.get_structure(atoms)


def collection(a: Any) -> Collection[Any]:
    """
    Wraps an object into a collection if the object is not a collection
    :param a: (Any) the object to wrap or a collection
    :return: (list of Any) the wrapped object
    """
    return [a] if isinstance(a, str) else (a if isinstance(a, (set, list, tuple)) else [a])


def unpack_single(a):
    return a[0] if len(a) == 1 else a


def add_method(cls: type) -> Callable:

    def decorator(func: Callable) -> Callable:

        @wraps(func)
        def wrapper(self, *args, **kwargs):
            return func(*args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator
