import numpy as np
from wmaee.codes.vasp_inputs import generate_potcar, automatic_kpoints, write_inputs
from wmaee.codes.vasp_runner import run_vasp
from wmaee.codes.vasp_outputs import parse_output
from wmaee.core.io import working_directory

# from wmaee.core.config import is_pmg_avail
# print(is_pmg_avail())

lat = 4.1*np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])

from ase import Atoms
struct = Atoms(
    "Cu",
    cell=lat,
    scaled_positions=[[0, 0 , 0]],
)

# from pymatgen.core import Structure
# struct = Structure(
#     lat,
#     ["Cu"],
#     [[0, 0 , 0]],
# )

# print(struct)
# print(type(struct))

# generate_potcar(struct, potcar_dir='/home/david/work/VASP_potentials/GGA-PBE/PAW_5.4')
# potcar = generate_potcar(struct, xc='LDA', potcar_mapping={'Cu': 'Cu_pv'})

incar = dict(NSW=0, IBRION=-1, EDIFF=1e-5, LWAVE=False, LCHARG=False)
kpoints = automatic_kpoints(50)

with working_directory('fcc-Cu'):
    write_inputs(struct=struct, incar=incar, 
                kpoints=kpoints, 
                xc='GGA-PBE', potcar_mapping={'Cu': 'Cu_pv'})
    run_vasp()
    output = parse_output()
    
print(output.final_energy)

# output = parse_output(directory='~/work/2021_03_NiTi_hydrides/data/raw/23-08_NiTiH_MD/conv_test/ortho_12H/LGL_5_LGA_5_PMS_1000', parse_xdatcar=True)
# print(len(output['ionic_relaxation']))
