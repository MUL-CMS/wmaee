import os
import numpy as np
import shutil

from pymatgen.io.vasp.inputs import Potcar, PotcarSingle
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure


def write_automatic_kpoints_file(length: float):
    '''
    This function writes a KPOINTS file which generates
    a regular Gamma-centered k-point mesh.

    Parameters
        ----------
        length : float
            Defines the length (R_k) that determines the subdevisions N1, N2, and N3.

    Returns
        -------
        File named KPOINTS

    '''
    with open('KPOINTS', 'w') as output:
        output.write('Fully automatic mesh\n')
        output.write('0 # indicates automatic generation scheme\n')
        output.write('Auto # fully automatic\n')
        output.write(f" {length} # length (R_k)")


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


def generate_potcar(structure: Structure, potcar_dir: str, potcar_mapping: dict=None, write: bool=True) -> Potcar:
    """
    Generate a POTCAR file for VASP calculations based on the elements in the structure.

    Parameters:
        structure : pymatgen.Structure
            The structure for which POTCAR is to be generated.
        potcar_dir : str
            The directory where the VASP POTCAR files are located.
        potcar_mapping : dict (optional, default: None)
            An optional mapping of species to POTCARs, e.g. {'Fe': 'Fe_pv'}. 
            For elements not explicitly specifies, a default mappig {'X': 'X'}
            will be applied.
        write : bool (optional, default: True)
            Whether to write the POTCAR (into current directory).
        

    Returns:
        pymatgen.Potcar:
            A Potcar object representing the POTCAR files for the elements in the structure.
    """
    
    
    # get list of species as in in POSCAR
    sp = [structure.sites[0].species_string]
    for s in structure.sites[1:]:
        if s.species_string != sp[-1]:
            sp.append(s.species_string)
    
    # get dictionary of all single POTCARs we need
    for s in set(sp):
        if s not in potcar_mapping.keys():
            potcar_mapping[s] = s
    potcars = {
       X: PotcarSingle.from_file(os.path.join(potcar_dir, potcar_mapping[X], 'POTCAR')) for X in sp
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
        with open("POTCAR", "w") as f:
            f.write(str(potcar))


def write_incar(incar_in: dict, magnetism=False, incar_out='INCAR'):
    '''
    Function writes an INCAR file.

    Parameters
        ----------
        incar_in : dict
            Takes a dictionary of key value pairs specifing the INCAR file.
        magnetism : boolean
            Needs to be completed for magnetism = True!
            If magnetism is True, then three additional parameters will be written
            to INCAR file: MAGMOM, ISPIN, LORBIT. Default = False.
        incar_out : str
            The name of the final INCAR file. Default = INCAR

    Returns
        -------
        File named INCAR.
    '''
    with open(incar_out, 'w') as f:
        for key, value in incar_in.items():
            f.write(key.replace(':', '=') + '=' + str(value) + '\n')
