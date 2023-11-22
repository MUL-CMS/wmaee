from wmaee.codes.lammps.lammps_potential import get_models, get_potentials
from wmaee.codes.lammps.lammps_runner import run_lammps
from wmaee.codes.lammps.lammps_outputs import parse_logfile
from ase.io.lammpsdata import write_lammps_data as _write_lammps_data
from ase.io.lammpsrun import read_lammps_dump as _read_lammps_dump
from typing import Union, IO, Any, List
from ase import Atoms

def write_lammps_data(dumpfile: Union[str, IO], atoms: Atoms, specorder: Optional[List[str]] = None, 
                      reduce_cell: bool = False, force_skew: bool = False, write_image_flags: bool = False,
                      masses: bool = False, velocities: bool = False, units: str = 'metal', 
                      atom_style: str = 'atomic'):
    """
    Write atomic structure data to a LAMMPS data file. Method forwarded from ase.io.lammpsdata.write_lammps_data.

    Parameters
    ----------
    dumpfile : Union[str, IO]
        File to which the output will be written.
    atoms : Atoms
        Atoms to be written.
    specorder : Optional[List[str]], optional
        Chemical symbols in the order of LAMMPS atom types, by default None
    force_skew : bool, optional
        Force to write the cell as a triclinic box (see LAMMPS documentation),
        by default False
    reduce_cell : bool, optional
        Whether the cell shape is reduced or not, by default False
    write_image_flags : bool, default False
        If True, the image flags indicating in which images of the periodic
        simulation box the atoms are, are written.
    masses : bool, optional
        Whether the atomic masses are written or not, by default False
    velocities : bool, optional
        Whether the atomic velocities are written or not, by default False
    units : str, optional
        LAMMPS units (see LAMMPS documentation), by default 'metal'
    atom_style : {'atomic', 'charge', 'full'}, optional
        LAMMPS atom style (see LAMMPS documentation), by default 'atomic'
    """
    _write_lammps_data(dumpfile, atoms, specorder=specorder, reduce_cell=reduce_cell, force_skew=force_skew, 
                       write_image_flags=write_image_flags, masses=masses, velocities=velocities, units=units,
                       atom_style=atom_style)
    

def read_lammps_dump(dump_file: Union[str, IO], index: Any = None, **kwargs) -> Union[Atoms, List[Atoms]]:
    """
    Read Atoms object(s) from a LAMMPS dump file. Method forwarded from ase.io.lammpsrun.read_lammps_dump.

    Parameters
    ----------
    dump_file : Union[str, IO]
        Name of the file to read from or a file descriptor.

    index : Any, optional
        Specifies which configuration to read. 
        Default is None, which returns the last configuration.
        
        Examples
        --------
        index=0
            First configuration
        index=-2
            Second to last
        index=':'
            All configurations
        index='-3:'
            Three last configurations
        index='::2'
            Even configurations
        index='1::2'
            Odd configurations

    Returns
    -------
    Union[Atoms, List[Atoms]]
        An Atoms object or a list of Atoms objects depending on the number of configurations read.
    """
    return _read_lammps_dump(dump_file, index=index)

__all__ = [get_models, get_potentials, run_lammps, parse_logfile, write_lammps_data, read_lammps_dump]