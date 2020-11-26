
try:
    import wmaee.extensions.external.parse_vasp_volumetric as parser
except ImportError:
    import logging
    logging.getLogger(__name__).info("Failed to load wmaee.extensions.external.parse_vasp_volumetric. "
                                     "Using Default Volumetric.from_file()")
else:
    from pymatgen.io.vasp import Poscar, VolumetricData
    from pymatgen.core import Lattice, Structure
    from types import MethodType

    def parse_file(file_name):
        raw = parser.parse_volumetric(file_name)
        structure = Structure(Lattice(raw["matrix"]), raw["species"], coords=raw["fcoords"])

        return Poscar(structure), dict(total=raw["charge_density"]), {}

    def _monkey_patch_parse_file():
        VolumetricData.parse_file = staticmethod(parse_file)

    _monkey_patch_parse_file()