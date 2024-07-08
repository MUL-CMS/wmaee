"""
Collection of helper functions for handling ML functionality of VASP.
"""

def ML_ABN_concat(files, output='ML_AB', overwrite=False):
    """Concatenates a list of `ML_AB` or `ML_ABN` files. Currently assumes that all files contain same chemistry!
    @TODO: check if the above is the case, and allow for merging various chemistries.

    Args:
        files (_type_): _description_
        output (str, optional): _description_. Defaults to 'ML_AB'.
        overwrite (bool, optional): _description_. Defaults to False.

    Raises:
        FileExistsError: _description_
    """
    
    from os.path import exists
    from os import remove, rename

    if exists(output) and not overwrite:
        raise FileExistsError(f'Output file {output} exists in the current folder and overwrite is not permitted.')
    else:
        out = open(output, 'w')

    s = 0
    for i, f in enumerate(files):
        # print(f'reading structures from {f}')
        with open(f, 'r') as src:
            # deal with header: first file copy header, other files ignore
            l = src.readline();
            while l.strip().split()[0] != 'Configuration':
                if i == 0:
                    out.writelines(l)
                l = src.readline()
            # done with headers, now go configuration by configuration
            if i > 0:
                out.writelines("**************************************************\n")
            while l:
                s += 1
                out.writelines(f"     Configuration num.{s:7d}\n")
                print(f"     Configuration num.{s:7d}")
                l = src.readline()
                while l and l.strip().split()[0] != 'Configuration':
                    out.writelines(l)
                    l = src.readline()
    out.close()

    # correct the total number of structure
    with open(output+'.tmp', 'w') as out:
        with open(output, 'r') as src:
            i = 1
            l = src.readline()
            while l:
                if i == 5:
                    out.writelines(f"{s:11d}\n")
                    print(f"{s:11d}")
                else:
                    out.writelines(l)
                l = src.readline()
                i += 1
    remove(output)
    rename(output+'.tmp', output)


def generate_ML_AB(input='OUTCAR', input_type='OUTCAR', output='ML_AB', overwrite=False, verbose=True):
    # from pymatgen.io.vasp import Vasprun
    from os.path import exists
    from monty.re import regrep
    from os import remove, rename
    
    # vrun = Vasprun(filename=vasprun, parse_dos=False, parse_eigen=False, parse_potcar_file=False, parse_projected_eigen=False)
    
    if not exists(input):
        raise FileExistsError(f'Input file {input} doesn\'t exist in the current folder.')
    if exists(output) and not overwrite:
        raise FileExistsError(f'Output file {output} exists in the current folder and overwrite is not permitted.')
    else:
        out = open(output, 'w')
        
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
            patterns=dict(num_species='ions per type\s+=\s+'+''.join(['(\d+)\s']*len(sp))),
            terminate_on_match=True,
            postprocess=lambda x: int(x)
            )
        num_sp = num_sp['num_species'][0][0]
    if verbose:
        print('num species: ', num_sp)
    if input_type=='OUTCAR':        
        num_sp = regrep(
            input,
            patterns=dict(masses='POMASS\s+=\s+'+''.join(['(\d+.\d+)\s']*len(sp))),
            terminate_on_match=True,
            postprocess=lambda x: float(x)
            )
        
    if verbose:
        print('num species: ', num_sp)
        
    
    