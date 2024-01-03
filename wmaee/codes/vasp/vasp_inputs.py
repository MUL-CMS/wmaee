from wmaee.core.config import Config
from wmaee.core.requirements import test_pmg
from wmaee.core.io import working_directory
from typing import Dict, Any, Optional, Union, Tuple
from ase import Atoms
from ase.io import write
import os
from shutil import copyfileobj

if test_pmg():
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.io.vasp import Potcar, PotcarSingle, Incar, Kpoints, Poscar


def automatic_kpoints(length: float, write: Optional[bool] = False) -> None:
    """
    Writes a KPOINTS file for generating a regular Gamma-centered k-point mesh.

    Parameters
    ----------
    length : float
        The length (R_k) that determines the subdivisions N1, N2, and N3.
    write : bool, optional
        If True, writes the KPOINTS file; if False, returns the KPOINTS content as a string.
        Default is False.

    Returns
    -------
    None or str
        If write is True, writes the KPOINTS file and returns None.
        If write is False, returns the KPOINTS content as a string.
    """
    kpts = 'Fully automatic mesh\n'
    kpts += '0 # indicates automatic generation scheme\n'
    kpts += 'Auto # fully automatic\n'
    kpts += f'{length} # length (R_k)'

    if write:
        with open('KPOINTS', 'w') as output:
            output.write(kpts)
    else:
        return kpts


def regular_kpoints(kpoints: Tuple[int, int, int] = (1, 1, 1),
                    shift: Tuple[float, float, float] = (0.0, 0.0, 0.0),
                    Gamma: bool = True,
                    write: bool = False) -> None:
    """
    Writes a KPOINTS file for generating a regular k-point mesh.

    Parameters
    ----------
    kpoints : tuple of int, optional
        Subdivisions N_1, N_2, and N_3 along the reciprocal lattice vectors.
        Default is (1, 1, 1).
    shift : tuple of float, optional
        Shift of the mesh (s_1, s_2, s_3).
        Default is (0.0, 0.0, 0.0).
    Gamma : bool, optional
        If True, generates a Gamma-centered mesh; if False, uses Monkhorst-Pack scheme.
        Default is True.
    write : bool, optional
        If True, writes the KPOINTS file; if False, returns the KPOINTS content as a string.
        Default is False.

    Returns
    -------
    None or str
        If write is True, writes the KPOINTS file and returns None.
        If write is False, returns the KPOINTS content as a string.
    """
    kpts = 'Regular k-point mesh\n'
    kpts += '0 # indicates automatic number of k-points\n'
    if Gamma:
        kpts += 'Gamma # generate a Gamma centered mesh\n'
    else:
        kpts += 'Monkhorst # Monkhorst-Pack scheme\n'
    kpts += f'{kpoints[0]} {kpoints[1]} {kpoints[2]} # subdivisions N_1, N_2, and N_3 along the reciprocal lattice vectors\n'
    kpts += f'{shift[0]} {shift[1]} {shift[2]} # shift of the mesh (s_1, s_2, s_3)'

    if write:
        with open('KPOINTS', 'w') as output:
            output.write(kpts)
    else:
        return kpts



def generate_potcar(struct: Union[Atoms, Any], 
                    potcar_dir: Optional[str] = None,
                    xc: Optional[str] = None, 
                    potcar_mapping: Optional[Dict[str, str]] = {},
                    write: Optional[bool] = True) -> Optional[Any]:
    """
    Generate a POTCAR file for VASP calculations based on the elements in the structure.

    Parameters
    ----------
    struct : Union[ase.Atoms, pymatgen.Structure]
        The structure for which POTCAR is to be generated. Can be either a pymatgen Structure or 
        an ASE Atoms object.
    potcar_dir : str, optional
        The directory where the VASP POTCAR files are located.
    xc : str, optional
        The exchange-correlation functional. Defaults to value in defined in
        `wmaee.config.yaml`.
    potcar_mapping : dict, optional
        An optional mapping of species to POTCARs, e.g. {'Fe': 'Fe_pv'}. 
        For elements not explicitly specified, a default mapping {'X': 'X'}
        will be applied.
    write : bool, optional
        Whether to write the POTCAR (into the current directory).

    Returns
    -------
    pymatgen.Potcar or None
        A Potcar object representing the POTCAR files for the elements in the structure 
        if pymatgen is available. Otherwise, returns None.
    """
    
    if potcar_dir == None:
        # No POTCAR dir specified, let's get it from config
        cfg = Config()
        root = cfg.get('applications').get('vasp').get('potential_path')
        if xc == None:
            xc = cfg.get('applications').get('vasp').get('potentials').get('default')
        xc_path = cfg.get('applications').get('vasp').get('potentials').get(xc)
        potcar_dir = os.path.join(root, xc_path)
    
    pmg = test_pmg()
    
    if pmg and isinstance(struct, Structure):
        # structure is pymatgen.core.Structure, convert to ase.Atoms
        struct = AseAtomsAdaptor.get_atoms(struct)
    if not isinstance(struct, Atoms):
        raise TypeError('Unknown structure type')
    
    # get list of species ordered as in structure
    sp = [struct[0].symbol]
    for s in struct[1:]:
        if s.symbol != sp[-1]:
            sp += [s.symbol]
    for s in set(sp):
        if s not in potcar_mapping.keys():
            potcar_mapping[s] = s
    if pmg:
        # let's use pymatgen's engine to get dictionary of all 
        # single POTCARs we need
        potcars = {
            X: PotcarSingle.from_file(os.path.join(potcar_dir, potcar_mapping[X], 'POTCAR')) 
            for X in sp
            }

        # construct POTCAR
        potcar = Potcar()
        functionals = []
        for X in sp:
            potcar.append(potcars[X])
            functionals.append(potcars[X].functional)
        potcar.functional = functionals[0]

        if write:
            # Write the POTCAR to a file
            with open('POTCAR', 'w') as outfile:
                outfile.write(str(potcar))
                
        return potcar
    
    # pymatgen is not available; should we write POTCAR? If so, let's just 
    # concatenate corresponding files
    if write:
        with open('POTCAR', 'w') as outfile:
            for X in sp:
                with open(os.path.join(potcar_dir, potcar_mapping[X], 'POTCAR'), 'r') as infile:
                    copyfileobj(infile, outfile)

    return None


import deprecation
from wmaee import __version__
@deprecation.deprecated(deprecated_in="2.1", removed_in="3.0",
                        current_version=__version__,
                        details="Use native pyMatGen or ASE functions")
def change_scale(filename: str, output_name: str, initial_scale: str, final_scale: str):
    '''
    This function changes the scale value in a POSCAR file. 
    Usefull for performing eV-calculations.

    Parameters
        ----------
        filename : str
            Name of the initial POSCAR file.
        output_name : str
            Name of final POSCAR file.
        initial_scale : str
            The initial scale written in filename (= initial POSCAR file).
        final_scale : str
            The scaling parameter that one wants to replace the initial scaling parameter with.


    Returns
        -------
        Modified POSCAR file.

    '''
    working_dir = os.getcwd()
    input_file = os.path.join(working_dir, filename)
    if not os.path.exists(input_file):
        raise FileNotFoundError(f'File {input_file} not found.')

    with open(input_file, 'r') as input_file:
        contents = input_file.read()

    # Look for the string specifying the scaling parameter in POSCAR
    # This will replace only the first occurance
    modified_contents = contents.replace(initial_scale, final_scale, 1)

    # Write the modified contents to the ouput file
    with open(output_name, 'w') as output_file:
        output_file.write(modified_contents)


import deprecation
from wmaee import __version__
@deprecation.deprecated(deprecated_in="2.1", removed_in="3.0",
                        current_version=__version__,
                        details="Use generate_potcar instead")
def concat_potcar_files(poscar_in='POSCAR', potcar_out='POTCAR', potcar_path=os.getcwd()):
    '''
    Creates a POTCAR file from individual atomic POTCAR files according to 
    the elements listed in POSCAR.

    Parameters
        ----------
        poscar_in : str
            Name of the initial POSCAR file to read in. Default = 'POSCAR
        potcar_out : str
            Name of final merged POTCAR file. Default = 'POTCAR'
        potcar_path : str
            The path were the individual POTCAR files are located. Default = current workind dir


    Returns
        -------
        Merged POTCAR file.
    '''

    potcars = []

    # list all files and directories in path
    files = os.listdir(potcar_path)

    for file in files:
        if 'POTCAR_' in file:
            potcars.append(file)

    # grep the elements from POSCAR
    elements_list = []
    with open(poscar_in, 'r') as poscar:
        lines = poscar.readlines()
        elements = lines[5].split(' ')
    for element in elements:
        element = element.strip('\n')
        elements_list.append(element)

    # Remove empty strings
    elements_list = [x for x in elements_list if x.strip()]

    # write POTCAR file in the same order of elements in POSCAR
    with open(potcar_out, 'wb') as outfile:
        for element in elements_list:
            found = False
            for potcar in potcars:
                if element in potcar:
                    found = True
                    print(f"{element} found in {potcar}")
                    with open(potcar, 'rb') as infile:
                        shutil.copyfileobj(infile, outfile)
                    break
            if not found:
                print(f"{element} not found in any potcar")
                

def write_incar(incar: Dict[str, Any], filename: str = 'INCAR') -> None:
    """
    Function writes an INCAR file.

    Parameters
    ----------
    incar : Dict[str, Any]
        A dictionary of key-value pairs specifying the INCAR file.
    filename : str, optional
        The name of the final INCAR file. Default is 'INCAR'.

    Returns
    -------
    None
        Writes the INCAR file.
    """
    with open(filename, 'w') as f:
        for key, value in incar.items():
            f.write(f"{key} = {value}\n")

            
# Using Any in type definitions to make it run even 
# in the case pymatgen is not installed
def write_inputs(struct: Union[Atoms, Any],
                 incar: Union[dict, Any],
                 kpoints: Optional[Union[None, str, Any]] = None,
                 potcar: Optional[Union[None, Any]] = None,
                 xc: Optional[Union[str, None]] = None,
                 potcar_mapping: Optional[Dict[str, str]] = {},
                 directory: Optional[Union[None, str]] = None) -> None:
    """
    Write VASP input files.
    Warning: If destination already contains any of the VASP inputs, they 
    get overwritten!

    Parameters
    ----------
    struct : Union[ase.Atoms, pymatgen.core.Structure, pymatgen.io.vasp.Poscar]
        Atomic structure information.
    incar : Union[dict, pymatgen.io.vasp.Incar]
        INCAR file information.
    kpoints : Optional[Union[None, str, pymatgen.io.vasp.Kpoints]], optional
        KPOINTS file information (default is None). If equal to None, INCAR must 
        contain k-points specification via KSPACING tag/
    potcar : Optional[Union[None, pymatgen.io.vasp.Potcar]], optional
        POTCAR file information (default is None and autogeneration based on 
        struct is used).
    xc : Optional[Union[str, None]], optional
        Exchange-correlation functional information (default is None and the 
        value is taken from wmaee.conf.yaml file).
    potcar_mapping : Optional[Dict[str, str]], optional
        Mapping of elements to POTCAR files (default is {}).
    directory : Optional[Union[None, str]], optional
        Directory to write files (default is None, e.g. current directory).

    Returns
    -------
    None
    """
    
    pmg = test_pmg()
    if directory == None:
        directory = '.'
    directory = os.path.expanduser(directory)
    with working_directory(directory):
        # POSCAR
        if isinstance(struct, Atoms):
            write(filename='POSCAR', images=struct, format='vasp')
        elif pmg and isinstance(struct, Poscar):
            struct.write_file(filename='POSCAR')
            struct = struct.structure # for POTCAR generation
        elif pmg and isinstance(struct, Structure):
            struct.to(filename='POSCAR', fmt='poscar')
        else:
            raise TypeError('Unknown structure type')
            
        # INCAR
        if isinstance(incar, dict):
            write_incar(incar)
        elif pmg and isinstance(incar, Incar):
            incar.write_file(filename='INCAR')
        else:
            raise TypeError('Unknown incar type')
        
        # KPOINTS
        if kpoints == None and 'KSPACING' not in incar:
            raise ValueError('kpoints must be specified when KSPACING is not in incar')
        elif isinstance(kpoints, str):
            with open('KPOINTS', 'w') as output:
                output.write(kpoints)
        elif pmg and isinstance(kpoints, Kpoints):
            kpoints.write_file(filename='KPOINTS')
            
        # POTCAR
        if pmg and isinstance(potcar, Potcar):
            with open('POTCAR', 'w') as outfile:
                outfile.write(str(potcar))
        else:
            generate_potcar(struct=struct, xc=xc, 
                            potcar_mapping=potcar_mapping, write=True)