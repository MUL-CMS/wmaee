from wmaee.codes.lammps.lammps_potential import get_models, get_potentials
from wmaee.codes.lammps.lammps_runner import run_lammps
from wmaee.codes.lammps.lammps_outputs import parse_logfile
from ase.io.lammpsdata import write_lammps_data as _write_lammps_data
from ase.io.lammpsrun import read_lammps_dump as _read_lammps_dump
from typing import Union, IO, Any, List, Optional
from ase import Atoms

def write_lammps_data(dumpfile: Union[str, IO], atoms: Atoms, specorder: Optional[List[str]] = None, 
                      force_skew: bool = False, velocities: bool = False, units: str = 'metal', 
                      atom_style: str = 'atomic', **kwargs):
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
    units : str, optional
        LAMMPS units (see LAMMPS documentation), by default 'metal'
    atom_style : {'atomic', 'charge', 'full'}, optional
        LAMMPS atom style (see LAMMPS documentation), by default 'atomic'
    """
    _write_lammps_data(dumpfile, atoms, specorder=specorder, force_skew=force_skew,
                       velocities=velocities, units=units, atom_style=atom_style, **kwargs)
    

def read_lammps_dump(dump_file: Union[str, IO], index: Any = slice(None), **kwargs) -> Union[Atoms, List[Atoms]]:
    """
    Read Atoms object(s) from a LAMMPS dump file. Method forwarded from ase.io.lammpsrun.read_lammps_dump.

    Parameters
    ----------
    dump_file : Union[str, IO]
        Name of the file to read from or a file descriptor.

    index : Any, optional
        Specifies which configuration to read. 
        Default is slice(None), which returns all configurations.
        
        Examples
        --------
        index=0
            First configuration
        index=-2
            Second to last
        index=slice(None)
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