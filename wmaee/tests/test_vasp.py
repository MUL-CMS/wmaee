import numpy as np
import wmaee.codes.vasp as vasp
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

# vasp.generate_potcar(struct, potcar_dir='/home/david/work/VASP_potentials/GGA-PBE/PAW_5.4')
# potcar = vasp.generate_potcar(struct, xc='LDA', potcar_mapping={'Cu': 'Cu_pv'})

incar = dict(NSW=0, IBRION=-1, EDIFF=1e-5, LWAVE=False, LCHARG=False)
kpoints = vasp.automatic_kpoints(50)

with working_directory('fcc-Cu'):
    vasp.write_inputs(struct=struct, incar=incar, 
                kpoints=kpoints, 
                xc='GGA-PBE', potcar_mapping={'Cu': 'Cu_pv'})
    vasp.run_vasp()
    output = vasp.parse_output()
    
print(output.final_energy)