import os
import numpy as np
from typing import Optional, Dict, Any, Union

from wmaee.core.data_structs import DotDict
from wmaee.core.requirements import test_pmg

if test_pmg():
    from pymatgen.io.vasp.outputs import Vasprun, Oszicar, Outcar, Xdatcar
    from pymatgen.io.vasp.inputs import Poscar

# for fallback solutions if pymatgen is not available
from ase.io.vasp import read_vasp, read_vasp_xml, read_vasp_out, read_vasp_xdatcar



# if possible, always use pymatgen for parsing
def parse_output(directory: Optional[str] = None,
                 POSCAR: str = 'POSCAR',
                 CONTCAR: str = 'CONTCAR',
                 OSZICAR: str = 'OSZICAR',
                 OUTCAR: str = 'OUTCAR',
                 vasprun: str = 'vasprun.xml',
                 XDATCAR: str = 'XDATCAR',                 
                 ase_atoms: bool = False,
                 parse_oszicar: bool = False,
                 parse_outcar: bool = False,
                 parse_vasprun: bool = True,
                 parse_vasprun_dos: bool = False,
                 parse_xdatcar: bool = False,
                 return_DocDict: bool = True,
                ) -> Union[DotDict, Dict[str, Any]]:
    """
    Parse VASP output files and return relevant information.

    Parameters
    ----------
    directory : str or None, optional
        The path to the directory containing VASP output files. If None, the current working directory is used.
    POSCAR : str, optional
        Name of the POSCAR file. Default is 'POSCAR'.
    CONTCAR : str, optional
        Name of the CONTCAR file. Default is 'CONTCAR'.
    OSZICAR : str, optional
        Name of the OSZICAR file. Default is 'OSZICAR'.
    OUTCAR : str, optional
        Name of the OUTCAR file. Default is 'OUTCAR'.
    vasprun : str, optional
        Name of the vasprun.xml file. Default is 'vasprun.xml'.
    XDATCAR : str, optional
        Name of the XDATCAR file. Default is 'XDATCAR'.
    ase_atoms : bool, optional
        Use ASE atoms instead of pymatgen structure if True. Default is False.
    parse_oszicar : bool, optional
        Parse ionic step energies from OSZICAR if True (only for pymatgen). Default is False.
    parse_outcar : bool, optional
        Parse final energies (and also magnetization if pymatgen is available) from OUTCAR if True. Default is False.
    parse_vasprun : bool, optional
        Parse vasprun.xml for ionic relaxation energies if True. Default is True.
    parse_vasprun_dos : bool, optional
        Parse DOS information from vasprun.xml if True (only if pymatgen is available). Default is False.
    parse_xdatcar : bool, optional
        Parse ionic relaxation steps from XDATCAR if True. Default is False.
    return_DocDict : bool, optional
        Return DotDict instead of a regular dictionary if True. Default is True.

    Returns
    -------
    DotDict or dict
        A DotDict or dictionary containing parsed information.
    """
    
    pmg = test_pmg()
    output_data = {}

    # Use the current working directory if vasp_dir is not provided
    if directory is None:
        directory = os.getcwd()
    directory = os.path.expanduser(directory)

    # Check if the directory exists
    if not os.path.isdir(directory):
        raise ValueError('Directory '+directory+' does not exist.')

    # Parse initial and final structures from POSCAR and CONTCAR
    filename = os.path.join(directory, POSCAR)
    if os.path.isfile(filename):
        if pmg and not ase_atoms:
            output_data['initial_structure'] = Poscar.from_file(filename=filename, check_for_POTCAR=False).structure
        else:
            output_data['initial_structure'] = read_vasp(filename)

    filename = os.path.join(directory, CONTCAR)
    if os.path.isfile(filename):
        if pmg and not ase_atoms:
            output_data['final_structure'] = Poscar.from_file(filename=filename, check_for_POTCAR=False).structure
        else:
            output_data['final_structure'] = read_vasp(filename)
            
    # Parse energies from OSZICAR if available - only pymatgen! - fast for energies
    oszicar_file = os.path.join(directory, OSZICAR)
    if parse_oszicar and pmg and os.path.isfile(oszicar_file):
        oszicar = Oszicar(oszicar_file)
        output_data['final_energy'] = oszicar.final_energy
        output_data['ionic_step_energies'] = [s['E0'] for s in oszicar.ionic_steps]
        
    # Parse data from OUTCAR
    if parse_outcar:
        filename = os.path.join(directory, OUTCAR)
        if os.path.isfile(filename):
            if pmg:
                # pymatgen parser
                outcar = Outcar(filename)
                output_data['total_magnetization'] = outcar.total_mag
                output_data['magnetization'] = outcar.magnetization
                output_data['final_energy'] = outcar.final_energy
            else:
                # ase (limited) parser
                outcar = read_vasp_out(filename, index=slice(None))
                energies = []
                forces = []
                stress = []
                for s in outcar:
                    energies.append(s.calc.get_potential_energy())
                    forces.append(np.array(s.calc.get_forces()))
                    stress.append(np.array(s.calc.get_stress()))
                output_data['final_energy'] = energies[-1]
                output_data['ionic_step_energies'] = energies
                output_data['ionic_step_forces'] = forces
                output_data['ionic_step_stress'] = stress
                output_data['final_stress'] = output_data['ionic_step_stress'][-1]
                output_data['final_forces'] = output_data['ionic_step_forces'][-1]
                # output_data['final_energy'] = outcar.calc.get_potential_energy()

    # Parse vasprun.xml for electronic structure information and DOS if enabled
    filename = os.path.join(directory, vasprun)
    if parse_vasprun and os.path.isfile(filename):
        if pmg: # pymatgen parser
            if parse_vasprun_dos:
                vrun = Vasprun(filename=filename, parse_eigen=False, parse_potcar_file=False, parse_projected_eigen=False, parse_dos=True)
                dos = vrun.complete_dos
                output_data['total_dos'] = dos
            else:
                vrun = Vasprun(filename=filename, parse_eigen=False, parse_potcar_file=False, parse_projected_eigen=False, parse_dos=False)
            output_data['final_energy'] = vrun.final_energy # first try to read final energy
            output_data['ionic_step_energies'] = [s['e_0_energy'] for s in vrun.ionic_steps]
            if 'forces' in vrun.ionic_steps[0].keys():
                output_data['ionic_step_forces'] = [np.array(s['forces']) for s in vrun.ionic_steps]
            if 'stress' in vrun.ionic_steps[0].keys():
                output_data['ionic_step_stress'] = [np.array(s['stress']) for s in vrun.ionic_steps]
        else: # ase parser
            vrun = read_vasp_xml(filename, index=slice(None))
            energies = []
            forces = []
            stress = []
            for s in vrun:
                energies.append(s.calc.get_potential_energy())
                forces.append(np.array(s.calc.get_forces()))
                stress.append(np.array(s.calc.get_stress()))
            output_data['final_energy'] = energies[-1]
            output_data['ionic_step_energies'] = energies
            output_data['ionic_step_forces'] = forces
            output_data['ionic_step_stress'] = stress
        output_data['final_stress'] = output_data['ionic_step_stress'][-1]
        output_data['final_forces'] = output_data['ionic_step_forces'][-1]
            
    # Parse ionic relaxation steps from XDATCAR
    filename = os.path.join(directory, XDATCAR)
    if parse_xdatcar and os.path.isfile(filename):
        if pmg:
            # pymatgen
            xdatcar = Xdatcar(filename)
            structures = xdatcar.structures
        else:
            # ase
            structures = read_vasp_xdatcar(filename, index=slice(None))
        output_data['ionic_relaxation'] = structures
            

    if return_DocDict:
        return DotDict(output_data)
    else:
        return output_data