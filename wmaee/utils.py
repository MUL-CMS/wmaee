from pymatgen.io.ase import AseAtomsAdaptor
from pyiron.atomistics.structure.atoms import ase_to_pyiron, pyiron_to_ase, ovito_to_pyiron, pyiron_to_ovito, \
    pyiron_to_pymatgen, pymatgen_to_pyiron, Atoms as IronAtoms
from typing import Any, Collection, Optional, Union
from pymatgen import Structure
from ase import Atoms as AseAtoms
from typing import Union, Callable, List
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


def add_method(cls: type, doc: Optional[Union[None, str]]) -> Callable:

    def decorator(func: Callable) -> Callable:

        @wraps(func)
        def wrapper(self, *args, **kwargs):
            return func(self, *args, **kwargs)
        if doc:
            wrapper.__doc__ = doc
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator


def animate(structures: List[Union[Structure, AseAtoms, IronAtoms]], spacefill: Optional[bool]=True, show_cell: Optional[bool]=True, stride: Optional[int]=1, center_of_mass: Optional[bool]=False, particle_size: Optional[float]=0.5):
    """
    Animates a structure and plays a movie using nglview
    :param structures: (list of Atoms or Structures) the images to visualize
    :param spacefill: (bool) wether to use spacefill style for the particles or not [default: True]
    :param show_cell: (bool) flag for displaying the unit cell [default: True]
    :param stride: (int) the stride, when you want to skip  images [default=1]
    :param center_of_mass: (bool) keep center of mass
    :param particle_size: (float) the default size of the particles
    :return: (nglview.view)
    """
    try:
        import nglview
    except ImportError:
        raise ImportError("The animate() function requires the package nglview to be installed")

    animation = nglview.show_asetraj([to_ase(s) for s in structures])
    if spacefill:
        animation.add_spacefill(radius_type='vdw', scale=0.5)
        animation.remove_ball_and_stick()
    else:
        animation.add_ball_and_stick()
    return animation
