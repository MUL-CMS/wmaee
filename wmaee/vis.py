from ase import Atoms as AseAtoms
from pyiron.atomistics.structure.atoms import Atoms as IronAtoms
from pymatgen import Structure
from wmaee.utils import to_pyiron, add_method as _hidden_add_method
from typing import Union, Optional, Any
from numpy import ndarray

try:
    from nglview import NGLWidget
except ImportError:
    raise


def view(structure: Union[Structure, AseAtoms, IronAtoms], show_cell: Optional[bool] = True,
         show_axes: Optional[bool] = True, camera: Optional[str] = 'orthographic', spacefill: Optional[bool] = True,
         particle_size: Optional[float] = 1.0, select_atoms: Optional[Union[None, ndarray]] = None,
         background: Optional[str] = 'white', color_scheme: Optional[Union[None, str]] = None,
         colors: Optional[Union[None, ndarray]] = None, scalar_field: Optional[Union[None, ndarray]] = None,
         scalar_start: Optional[Union[None, float]] = None, scalar_end: Optional[Union[None, float]] = None,
         scalar_cmap: Optional[Any] = None, vector_field: Optional[Union[None, ndarray]] = None,
         vector_color: Optional[Union[None, ndarray]] = None, custom_array: Optional[Union[None, ndarray]] = None,
         custom_3darray: Optional[Union[None, ndarray]] = None):
    """
    Plot3d relies on NGLView to visualize atomic structures. Here, we construct a string in the "protein database"
    ("pdb") format, then turn it into an NGLView "structure". PDB is a white-space sensitive format, so the
    string snippets are carefully formatted.
    The final widget is returned. If it is assigned to a variable, the visualization is suppressed until that
    variable is evaluated, and in the meantime more NGL operations can be applied to it to modify the visualization.

    :param show_cell: (bool) Whether or not to show the frame. (Default is True.)
    :param show_axes: (bool) Whether or not to show xyz axes. (Default is True.)
    :param camera: (str) 'perspective' or 'orthographic'. (Default is 'perspective'.)
    :param spacefill: (bool) Whether to use a space-filling or ball-and-stick representation. (Default is True, use space-filling atoms.)
    :param particle_size: (float) Size of the particles. (Default is 1.)
    :param select_atoms: (numpy.ndarray) Indices of atoms to show, either as integers or a boolean array mask. (Default is None, show all atoms.)
    :param background: (str) Background color. (Default is 'white'.)
    :param color_scheme: (str) NGLView color scheme to use. (Default is None, color by element.)
    :param colors: (numpy.ndarray) A per-atom array of HTML color names or hex color codes to use for atomic colors. (Default is None, use coloring scheme.)
    :param scalar_field: (numpy.ndarray) Color each atom according to the array value (Default is None, use coloring scheme.)
    :param scalar_start: (float) The scalar value to be mapped onto the low end of the color map (lower values are clipped). (Default is None, use the minimum value in `scalar_field`.)
    :param scalar_end: (float) The scalar value to be mapped onto the high end of the color map (higher values are clipped). (Default is None, use the maximum value in `scalar_field`.)
    :param scalar_cmap: (matplotlib.cm) The colormap to use. (Default is None, giving a blue-red divergent map.)
    :param vector_field: (numpy.ndarray) Add vectors (3 values) originating at each atom. (Default is None, no vectors.)
    :param vector_color: (numpy.ndarray) Colors for the vectors (only available with vector_field). (Default is None, vectors are colored by their direction.)
    Possible NGLView color schemes
          " ", "picking", "random", "uniform", "atomindex", "residueindex",
          "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
          "hydrophobicity", "value", "volume", "occupancy"
    :return (nglview.NGLWidget): The NGLView widget itself, which can be operated on further or viewed as-is.
    Warnings:
        * Many features only work with space-filling atoms (e.g. coloring by a scalar field).
        * The colour interpretation of some hex codes is weird, e.g. 'green'.
    """
    return to_pyiron(structure).plot3d(show_cell=show_cell, show_axes=show_axes, camera=camera, spacefill=spacefill,
                                       particle_size=particle_size, select_atoms=select_atoms, background=background,
                                       color_scheme=color_scheme, colors=colors, scalar_field=scalar_field,
                                       scalar_start=scalar_start, scalar_end=scalar_end, scalar_cmap=scalar_cmap,
                                       vector_color=vector_color, vector_field=vector_field,
                                       custom_3darray=custom_3darray, custom_array=custom_array)

_non_forward_view = view

@_hidden_add_method(Structure, view.__doc__)
def view(self, show_cell: Optional[bool] = True,
         show_axes: Optional[bool] = True, camera: Optional[str] = 'orthographic', spacefill: Optional[bool] = True,
         particle_size: Optional[float] = 1.0, select_atoms: Optional[Union[None, ndarray]] = None,
         background: Optional[str] = 'white', color_scheme: Optional[Union[None, str]] = None,
         colors: Optional[Union[None, ndarray]] = None, scalar_field: Optional[Union[None, ndarray]] = None,
         scalar_start: Optional[Union[None, float]] = None, scalar_end: Optional[Union[None, float]] = None,
         scalar_cmap: Optional[Any] = None, vector_field: Optional[Union[None, ndarray]] = None,
         vector_color: Optional[Union[None, ndarray]] = None, custom_array: Optional[Union[None, ndarray]] = None,
         custom_3darray: Optional[Union[None, ndarray]] = None):
    return to_pyiron(self).plot3d(show_cell=show_cell, show_axes=show_axes, camera=camera, spacefill=spacefill,
                                       particle_size=particle_size, select_atoms=select_atoms, background=background,
                                       color_scheme=color_scheme, colors=colors, scalar_field=scalar_field,
                                       scalar_start=scalar_start, scalar_end=scalar_end, scalar_cmap=scalar_cmap,
                                       vector_color=vector_color, vector_field=vector_field,
                                       custom_3darray=custom_3darray, custom_array=custom_array)

@_hidden_add_method(AseAtoms, view.__doc__)
def view(self, show_cell: Optional[bool] = True,
         show_axes: Optional[bool] = True, camera: Optional[str] = 'orthographic', spacefill: Optional[bool] = True,
         particle_size: Optional[float] = 1.0, select_atoms: Optional[Union[None, ndarray]] = None,
         background: Optional[str] = 'white', color_scheme: Optional[Union[None, str]] = None,
         colors: Optional[Union[None, ndarray]] = None, scalar_field: Optional[Union[None, ndarray]] = None,
         scalar_start: Optional[Union[None, float]] = None, scalar_end: Optional[Union[None, float]] = None,
         scalar_cmap: Optional[Any] = None, vector_field: Optional[Union[None, ndarray]] = None,
         vector_color: Optional[Union[None, ndarray]] = None, custom_array: Optional[Union[None, ndarray]] = None,
         custom_3darray: Optional[Union[None, ndarray]] = None):
    return to_pyiron(self).plot3d(show_cell=show_cell, show_axes=show_axes, camera=camera, spacefill=spacefill,
                                       particle_size=particle_size, select_atoms=select_atoms, background=background,
                                       color_scheme=color_scheme, colors=colors, scalar_field=scalar_field,
                                       scalar_start=scalar_start, scalar_end=scalar_end, scalar_cmap=scalar_cmap,
                                       vector_color=vector_color, vector_field=vector_field,
                                       custom_3darray=custom_3darray, custom_array=custom_array)

@_hidden_add_method(IronAtoms, view.__doc__)
def view(self, show_cell: Optional[bool] = True,
         show_axes: Optional[bool] = True, camera: Optional[str] = 'orthographic', spacefill: Optional[bool] = True,
         particle_size: Optional[float] = 1.0, select_atoms: Optional[Union[None, ndarray]] = None,
         background: Optional[str] = 'white', color_scheme: Optional[Union[None, str]] = None,
         colors: Optional[Union[None, ndarray]] = None, scalar_field: Optional[Union[None, ndarray]] = None,
         scalar_start: Optional[Union[None, float]] = None, scalar_end: Optional[Union[None, float]] = None,
         scalar_cmap: Optional[Any] = None, vector_field: Optional[Union[None, ndarray]] = None,
         vector_color: Optional[Union[None, ndarray]] = None, custom_array: Optional[Union[None, ndarray]] = None,
         custom_3darray: Optional[Union[None, ndarray]] = None):
    return to_pyiron(self).plot3d(show_cell=show_cell, show_axes=show_axes, camera=camera, spacefill=spacefill,
                                       particle_size=particle_size, select_atoms=select_atoms, background=background,
                                       color_scheme=color_scheme, colors=colors, scalar_field=scalar_field,
                                       scalar_start=scalar_start, scalar_end=scalar_end, scalar_cmap=scalar_cmap,
                                       vector_color=vector_color, vector_field=vector_field,
                                       custom_3darray=custom_3darray, custom_array=custom_array)