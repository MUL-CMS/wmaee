import numpy as np
from ase import Atoms
from pymatgen.core import Structure
from wmaee.codes.vasp_inputs import generate_potcar, automatic_kpoints, write_inputs
from wmaee.codes.vasp_runner import run_vasp
from wmaee.core.io import working_directory

lat = 4.1*np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])

struct = Atoms(
    "Cu",
    cell=lat,
    scaled_positions=[[0, 0 , 0]],
)
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