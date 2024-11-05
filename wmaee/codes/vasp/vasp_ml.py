"""
Collection of helper functions for handling ML functionality of VASP.
"""

from typing import Callable
from typing import List, Union

def ML_ABN_concat(files: List[str], output: str = 'ML_AB', overwrite: bool = False, verbose: bool = True) -> None:
    """
    Concatenates a list of `ML_AB` or `ML_ABN` files. This function merges headers and structures
    from multiple files into a single output file.

    Parameters
    ----------
    files : List[str]
        List of paths to input files that should be concatenated.
    output : str, optional
        Path for the output file, by default 'ML_AB'.
    overwrite : bool, optional
        Whether to overwrite the output file if it already exists, by default False.
    verbose : bool, optional
        If True, prints details of the merging process, by default True.

    Raises
    ------
    FileExistsError
        If the output file exists and overwrite is set to False.

    Notes
    -----
    This function currently assumes that all input files contain the same chemistry.
    """
    
    from os.path import exists
    from os import remove, rename

    if exists(output) and not overwrite:
        raise FileExistsError(f'Output file {output} exists in the current folder and overwrite is not permitted.')
    else:
        out = open(output+'.tmp_wmaee_tmp', 'w')

    s = 0
    sp = list()
    max_types = 0
    max_atoms = 0
    max_species = 0
    masses = dict()
    for i, f in enumerate(files):
        if verbose:
            print(f'reading structures from {f}')
        with open(f, 'r') as src:
            # deal with header: collect merged information
            l = src.readline();
            while l.strip().split()[0] != 'Configuration':                
                if 'The atom types' in l.strip():
                    l = src.readline();
                    l = src.readline();
                    sp += l.strip().split()
                    if verbose:
                        print(f'  found species: {sp}')
                if 'number of atom type' in l.strip():
                    l = src.readline();
                    l = src.readline();
                    max_at = int(l.strip())
                    if verbose:
                        print(f'  maximum number of atom typess: {max_at}')
                    max_types = max(max_types, max_at)                        
                if 'atoms per system' in l.strip():
                    l = src.readline();
                    l = src.readline();
                    max_at = int(l.strip())
                    if verbose:
                        print(f'  maximum number of atoms: {max_at}')
                    max_atoms = max(max_atoms, max_at)
                if 'atoms per atom type' in l.strip():
                    l = src.readline();
                    l = src.readline();
                    max_at = int(l.strip())
                    if verbose:
                        print(f'  maximum number of atoms per type: {max_at}')
                    max_species = max(max_species, max_at)
                if 'Atomic mass' in l.strip():
                    l = src.readline();
                    m = ''
                    l = src.readline();
                    while '*********************' not in l.strip():
                        m += l
                        l = src.readline();
                    for X, M in zip(sp, m.strip().split()):
                        masses[X] = float(M)
                        if verbose:
                            print(f'  mass of {X}: {float(M)} a.u.')              
                l = src.readline();
            # done with headers, now go configuration by configuration            
            out.writelines("**************************************************\n")
            while l:
                s += 1
                out.writelines(f"     Configuration num.{s:7d}\n")
                if verbose:
                    print(f"     Configuration num.{s:7d}")
                l = src.readline()
                while l and l.strip().split()[0] != 'Configuration':
                    out.writelines(l)
                    l = src.readline()
    out.close()
    
    sp = list(set(sp))

    # write header    
    with open(output, 'w') as out:        
        _SEP_star  = "**************************************************\n"
        _SEP_minus = "--------------------------------------------------\n"
        _SEP_equal = "==================================================\n"
        out.writelines(' 1.0 Version\n'+_SEP_star)
        out.writelines('     The number of configurations\n'+_SEP_minus+f'{s:11d}\n'+_SEP_star)
        out.writelines('     The maximum number of atom type\n'+_SEP_minus+f'{max_types:8d}\n'+_SEP_star)
        out.writelines('     The atom types in the data file\n'+_SEP_minus+'     ')
        for i, el in enumerate(sp):
            out.writelines(el.ljust(3))
            if i%3 == 2:
                out.writelines('\n     ')
        out.writelines('\n'+_SEP_star)
        out.writelines('     The maximum number of atoms per system\n'+_SEP_minus+f'{max_atoms:8d}\n'+_SEP_star)
        out.writelines('     The maximum number of atoms per atom type\n'+_SEP_minus+f'{max_species:8d}\n'+_SEP_star)
        out.writelines('     Reference atomic energy (eV)\n'+_SEP_minus)
        for i, _ in enumerate(sp):
            out.writelines('   '+f'{0:18.16f}'+'     ') # dummy energy
            if i%3 == 2:
                out.writelines('\n')
        out.writelines('\n'+_SEP_star)
        out.writelines('     Atomic mass\n'+_SEP_minus)        
        for i, X in enumerate(sp):
            out.writelines('   '+f'{masses[X]:18.16f}'+'     ')
            if i%3 == 2:
                out.writelines('\n')
        out.writelines('\n'+_SEP_star)
        out.writelines('     The numbers of basis sets per atom type\n'+_SEP_minus+'    ')
        for i, _ in enumerate(sp):
            out.writelines('     1')
            if i%3 == 2:
                out.writelines('\n')
        out.writelines('\n')
        for i, el in enumerate(sp):
            out.writelines(_SEP_star+f'     Basis set for {el}\n'+_SEP_minus+'          1      1\n')
        
        with open(output+'.tmp_wmaee_tmp', 'r') as src:
            l = src.readline()
            while l:
                out.writelines(l)
                l = src.readline()
              
    remove(output+'.tmp_wmaee_tmp')
    
    
    
def _read_table_pattern(filename: str, header_pattern: str, row_pattern: str,
        footer_pattern: str, postprocess: Callable = str, first_one_only: bool = False,
    ) -> list:
    """
    Parse table-like data from a file. A table comprises three parts: header,
    main body, and footer. All the data matching the "row pattern" in the main body
    will be returned.

    Parameters
    ----------
    filename : str
        The name of the file to read the table from.
    header_pattern : str
        The regular expression pattern that matches the table header.
        This pattern should match all the text immediately before the main body of the table.
        For multiple sections table, match the text until the section of interest.
        MULTILINE and DOTALL options are enforced, so the "." meta-character will also match "\n" in this section.
    row_pattern : str
        The regular expression that matches a single line in the table. Capture interested fields using regular expression groups.
    footer_pattern : str
        The regular expression that matches the end of the table, e.g., a long dash line.
    postprocess : Callable, optional
        A post-processing function to convert all matches. Defaults to str, i.e., no change.
    first_one_only : bool, optional
        Only the first occurrence of the table will be parsed and the parsing procedure will stop. 
        The enclosing list will be removed, i.e., only a single table will be returned. Incompatible with last_one_only.

    Returns
    -------
    List[Union[List, Dict]]
        A list of tables. A table is a list of rows. A row is either a list of attribute values (if capturing groups are defined without names in row_pattern), 
        or a dictionary (if named capturing groups are defined by row_pattern).

    Notes
    -----
    This function is adapted from pymatgen.io.vasp.outputs.py.
    """
    
    from monty.io import zopen
    import re
    
    with zopen(filename, mode='rt') as file:
        text = file.read()
    table_pattern_text = header_pattern + r"\s*^(?P<table_body>(?:\s+" + row_pattern + r")+)\s+" + footer_pattern
    table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
    rp = re.compile(row_pattern)
    tables: list[list] = []
    for mt in table_pattern.finditer(text):
        table_body_text = mt.group("table_body")
        table_contents = []
        for line in table_body_text.split("\n"):
            ml = rp.search(line)
            # Skip empty lines
            if not ml:
                continue
            d = ml.groupdict()
            if len(d) > 0:
                processed_line: dict | list = {k: postprocess(v) for k, v in d.items()}
            else:
                processed_line = [postprocess(v) for v in ml.groups()]
            table_contents.append(processed_line)
        tables.append(table_contents)
        if first_one_only:
            break
    return tables
    
    

def generate_ML_AB(input: str = 'OUTCAR', input_type: str = 'OUTCAR', output: str = 'ML_AB', 
                   overwrite: bool = False, verbose: bool = True, system_name: str = 'parsed AIMD') -> None:   
    """
    Generate ML_AB formatted output from VASP calculation outputs.

    Parameters
    ----------
    input : str, optional
        Path to the input file, by default 'OUTCAR'
    input_type : str, optional
        Type of the input file, either 'OUTCAR' or 'vasprun', by default 'OUTCAR'
    output : str, optional
        Path to the output file, by default 'ML_AB'
    overwrite : bool, optional
        Flag to overwrite the existing output file, by default False
    verbose : bool, optional
        Flag to print detailed output, by default True
    system_name : str, optional
        Name of the system, by default 'parsed AIMD'

    Raises
    ------
    Exception
        If input_type is not 'OUTCAR' or 'vasprun'
    FileExistsError
        If the input file does not exist
    """
    
    from os.path import exists
    from monty.re import regrep
    import numpy as np
    
     # Validate input_type
    if input_type not in ['OUTCAR', 'vasprun']:
        raise Exception(f'Unknown input file type, allowed values are `OUTCAR` or `vasprun`')

    # Check if input file exists
    if not exists(input):
        raise FileExistsError(f'Input file {input} doesn\'t exist in the current folder.')

    # Read species, numbers, and masses from the input file
    if input_type=='vasprun':
        chemistry = _read_table_pattern(
            filename=input,            
            header_pattern=r"\s*<field type=\"string\">pseudopotential<\/field>\s*<set>",
            row_pattern=r"\s*<rc><c>\s*(\d*)<\/c><c>\s*(\w*)\s*<\/c><c>\s*(\d*.\d*)<\/c><c>.*",
            footer_pattern=r"\s*<\/set>",
            first_one_only=True
        )[0]
        sp = list()
        num_sp = list()
        masses = list()
        for n, s, m in chemistry:
            sp.append(s)
            num_sp.append(int(n))
            masses.append(float(m))
    if input_type=='OUTCAR':
        sp = regrep(
            input, 
            patterns=dict(species='TITEL\s+=\s\w+\s(\w+)')
            )
        sp = [x[0][0] for x in sp['species']]
    if verbose:
        print('found species: ', sp)
    if input_type=='OUTCAR':
        num_sp = regrep(
            input,
            patterns=dict(num_species='ions per type\s+=\s+'+''.join(['(\d+)\s*']*len(sp))),
            terminate_on_match=True,
            postprocess=lambda x: int(x)
            )
        # print('ions per type\s+=\s+'+''.join(['(\d+)\s*']*len(sp)))
        num_sp = num_sp['num_species'][0][0]
    if verbose:
        print('num species: ', num_sp)
    if input_type=='OUTCAR':
        masses = regrep(
            input,
            patterns=dict(masses='POMASS\s+=\s+'+''.join(['(\d+.\d+)\s*']*len(sp))),
            terminate_on_match=True,
            postprocess=lambda x: float(x)
            )
        masses = masses['masses'][0][0]
    if verbose:
        print('masses: ', masses)
        
    # Read lattice vectors from the input file
    if input_type == 'OUTCAR':
        lattices = _read_table_pattern(
            filename=input,
            header_pattern=r"\sdirect lattice vectors\s+reciprocal lattice vectors",
            row_pattern=r"\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+",
            footer_pattern=r"\n\s+length of vectors",
            postprocess=lambda x: float(x)
        )[2:]  # Skip initial lattice
    else:
        lattices = _read_table_pattern(
            filename=input,
            header_pattern=r"\s*<varray name=\"basis\"\s*>",
            row_pattern=r"\s*<v>\s*([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s<\/v>\s*",
            footer_pattern=r"\s*<\/varray>",
            postprocess=lambda x: float(x)
        )[2:]  # Skip initial lattice
    if len(lattices[-1]) != 3:
        lattices = lattices[:-1]  # Remove invalid lattice
    if verbose:
        print('found lattices: ', len(lattices))
        
    # Read energies and stresses from the input fi
    if input_type=='OUTCAR':
        energies_stresses = regrep(
            input,
            patterns=dict(
                energies='free  energy\s+TOTEN\s+=\s+([+-]?\d+.\d+)\seV',
                stresses='in kB\s+([+-]?\d+.\d+)\s+([+-]?\d+.\d+)\s+([+-]?\d+.\d+)\s+([+-]?\d+.\d+)\s+([+-]?\d+.\d+)\s+([+-]?\d+.\d+)'
                ),
            terminate_on_match=False,
            postprocess=lambda x: float(x)
            )
        energies = [e[0][0] for e in energies_stresses['energies']]
        stresses = [s[0] for s in energies_stresses['stresses']]
    else:
        energies = _read_table_pattern(
            filename=input,
            header_pattern=r"\s*<energy>",
            row_pattern=r"\s*<i name=\"\w*\">\s*([+-]?\d+\.\d+)\s+<\/i>\s*",
            footer_pattern=r"\s*<i name=\"kinetic\">",
            postprocess=lambda x: float(x))        
        for i, e in enumerate(energies):
            energies[i] = e[0][0]
        stresses = _read_table_pattern(
            filename=input,
            header_pattern=r"\s*<varray name=\"stress\"\s*>",
            row_pattern=r"\s*<v>\s*([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s<\/v>\s*",
            footer_pattern=r"\s*<\/varray>",
            postprocess=lambda x: float(x))
        for i, s in enumerate(stresses):
            stresses[i] = [s[0][0], s[1][1], s[2][2], s[0][1], s[0][2], s[1][2]]
    if verbose:
        print('found energies: ', len(energies))
        print('found stresses: ', len(stresses))
        
    # Read positions and forces from the input file
    if input_type=='OUTCAR':
        positions_and_forces = _read_table_pattern(
            filename=input,
            header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
            row_pattern=r"\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
            footer_pattern=r"\s--+",
            postprocess=lambda x: float(x))
        if len(positions_and_forces[-1]) != sum(num_sp):
            # invalid number of positions and forces of the last structure, remove it
            positions_and_forces = positions_and_forces[:-1]
        positions = [[x[:3] for x in ionic_step] for ionic_step in positions_and_forces]
        forces = [[x[3:] for x in ionic_step] for ionic_step in positions_and_forces]
    else:
        positions = _read_table_pattern(
            filename=input,
            header_pattern=r"\s*<varray name=\"positions\"\s*>",
            row_pattern=r"\s*<v>\s*([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s<\/v>\s*",
            footer_pattern=r"\s*<\/varray>",
            postprocess=lambda x: float(x))[2:] # the first lattices are printed during initialization, not after ionic steps
        forces = _read_table_pattern(
            filename=input,
            header_pattern=r"\s*<varray name=\"forces\"\s*>",
            row_pattern=r"\s*<v>\s*([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s<\/v>\s*",
            footer_pattern=r"\s*<\/varray>",
            postprocess=lambda x: float(x))
    if verbose:
        print('found positions: ', len(positions))
        print('found forces: ', len(forces))
    
    # Take the minimum number of data points across all lists    
    take_only = min(len(lattices), len(positions), len(energies), len(stresses))
    lattices = lattices[:take_only]
    energies = energies[:take_only]
    positions = positions[:take_only]
    
    # Write ML_AB file
    out = open(output, 'w')
    _SEP_star  = "**************************************************\n"
    _SEP_minus = "--------------------------------------------------\n"
    _SEP_equal = "==================================================\n"
    out.writelines(' 1.0 Version\n'+_SEP_star)
    out.writelines('     The number of configurations\n'+_SEP_minus+f'{len(lattices):11d}\n'+_SEP_star)
    out.writelines('     The maximum number of atom type\n'+_SEP_minus+f'{len(num_sp):8d}\n'+_SEP_star)
    out.writelines('     The atom types in the data file\n'+_SEP_minus+'     ')
    for i, el in enumerate(sp):
        out.writelines(el.ljust(3))
        if i%3 == 2:
            out.writelines('\n     ')
    out.writelines('\n'+_SEP_star)
    out.writelines('     The maximum number of atoms per system\n'+_SEP_minus+f'{sum(num_sp):8d}\n'+_SEP_star)
    out.writelines('     The maximum number of atoms per atom type\n'+_SEP_minus+f'{max(num_sp):8d}\n'+_SEP_star)
    out.writelines('     Reference atomic energy (eV)\n'+_SEP_minus)
    for i, _ in enumerate(sp):
        out.writelines('   '+f'{0:18.16f}'+'     ') # dummy energy
        if i%3 == 2:
            out.writelines('\n')
    out.writelines('\n'+_SEP_star)
    out.writelines('     Atomic mass\n'+_SEP_minus)
    for i, m in enumerate(masses):
        out.writelines('   '+f'{m:18.16f}'+'     ')
        if i%3 == 2:
            out.writelines('\n')
    out.writelines('\n'+_SEP_star)
    out.writelines('     The numbers of basis sets per atom type\n'+_SEP_minus+'    ')
    for i, _ in enumerate(sp):
        out.writelines('     1')
        if i%3 == 2:
            out.writelines('\n')
    out.writelines('\n')
    for i, el in enumerate(sp):
        out.writelines(_SEP_star+f'     Basis set for {el}\n'+_SEP_minus+'          1      1\n')
    for i in range(len(lattices)):
        out.writelines(_SEP_star)
        out.writelines(f'     Configuration num.{i+1:7d}\n'+_SEP_equal)
        out.writelines('     System name\n'+_SEP_minus+'     '+system_name+'\n'+_SEP_equal)
        out.writelines('     The number of atom types\n'+_SEP_minus+f'{len(num_sp):8d}\n'+_SEP_equal)
        out.writelines('     The number of atoms\n'+_SEP_minus+f'{sum(num_sp):8d}\n'+_SEP_star)
        out.writelines('     Atom types and atom numbers\n'+_SEP_minus)
        for el, n in zip(sp, num_sp):
            out.writelines('     '+el.ljust(3)+f'{n:6d}\n')
        out.writelines(_SEP_equal)
        out.writelines('     Primitive lattice vectors (ang.)\n'+_SEP_minus)
        lattice = lattices[i]
        for l in lattice:
            out.writelines('  '+'  '.join([f'{x:18.16f}' for x in l])+'\n')
        out.writelines(_SEP_equal)
        out.writelines('     Atomic positions (ang.)\n'+_SEP_minus)
        if input_type == 'OUTCAR':
            for l in positions[i]:
                out.writelines('  '+'  '.join([f'{x:18.16f}' for x in l])+'\n')
        else:
            # need to convert to Cartesian coordinates for vasprun.xml inputs
            for l in positions[i]:
                out.writelines('  '+'  '.join([f'{x:18.16f}' for x in np.matmul(l, np.array(lattice))])+'\n')
        out.writelines(_SEP_equal)
        out.writelines('     Total energy (eV)\n'+_SEP_minus+f'  {energies[i]}\n'+_SEP_equal)
        out.writelines('     Forces (eV ang.^-1)\n'+_SEP_minus)
        for l in forces[i]:
            out.writelines('  '+'  '.join([f'{x:18.16f}' for x in l])+'\n')
        out.writelines(_SEP_equal)
        out.writelines('     Stress (kbar)\n'+_SEP_minus+'     XX YY ZZ\n'+_SEP_minus)
        out.writelines('  '+'  '.join([f'{x:18.16f}' for x in stresses[i][:3]])+'\n'+_SEP_minus)
        out.writelines('     XY YZ ZX\n'+_SEP_minus)
        out.writelines('  '+'  '.join([f'{x:18.16f}' for x in stresses[i][3:]])+'\n')
