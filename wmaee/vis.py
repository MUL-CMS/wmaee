from ase import Atoms
from pymatgen import Structure
from wmaee.utils import pymatgen_to_ase

def view(structure, spacefill=True, show_cell=True, camera='perspective', particle_size=0.5, background='white', color_scheme='element', show_axes=True):
    """
    Constructs a nglview view to display a structure
    :param structure: (pymatgen.Structure, ase.Atoms) the structure to display
    :param spacefill: (bool) to set the atoms size to spacefilling (default: True)
    :param show_cell:  (bool) wether to draw the unit cell or not (default: True)
    :param camera: (str) which camera projections to use 'perspective' or 'orthographic' (default: 'perspective')
    :param particle_size: (float) the size of the atoms (default: 0.5)
    :param background: (str) the name of the background color (default: 'white')
    :param color_scheme: (str) the name of the coloring scheme. Please refer to nglview documentation (default: 'element')
    :param show_axes: (bool) wether to draw the coordinate system of the unit cell (default: True)
    :return: (nglview.View) view wrapper
    """
    try:
        import nglview
    except ImportError:
        raise ImportError('nglview is needed')
    if isinstance(structure, Atoms):
        atoms = structure
    elif isinstance(structure, Structure):
        atoms = pymatgen_to_ase(structure)
    else:
        raise TypeError
    view_  = nglview.show_ase(atoms)
    if spacefill:
        view_.add_spacefill(radius_type='vdw', color_scheme=color_scheme, radius=particle_size)
        # view.add_spacefill(radius=1.0)
        view_.remove_ball_and_stick()
    else:
        view_.add_ball_and_stick()
    if show_cell:
        if atoms.cell is not None:
            view_.add_unitcell()
    if show_axes:
        view_.shape.add_arrow([-2, -2, -2], [2, -2, -2], [1, 0, 0], 0.5)
        view_.shape.add_arrow([-2, -2, -2], [-2, 2, -2], [0, 1, 0], 0.5)
        view_.shape.add_arrow([-2, -2, -2], [-2, -2, 2], [0, 0, 1], 0.5)
    if camera != 'perspective' and camera != 'orthographic':
        print('Only perspective or orthographic is permitted')
        return None
    view_.camera = camera
    view_.background = background
    return view_
