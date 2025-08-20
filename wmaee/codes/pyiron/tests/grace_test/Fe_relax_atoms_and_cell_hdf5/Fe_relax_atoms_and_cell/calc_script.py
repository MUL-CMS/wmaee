from tensorpotential.calculator import grace_fm
from ase.io import read
from ase.io.cif import write_cif
from ase.io.trajectory import Trajectory
from ase.optimize import FIRE, BFGS
from ase.constraints import ExpCellFilter as ECF
initial_structure = read("structure.cif")
calc = grace_fm("GRACE-2L-OMAT")
initial_structure.calc = calc
relaxation = FIRE(ECF(atoms=initial_structure, ), trajectory="relax.traj").run(fmax=0.1, steps=500)
relaxed_structure = initial_structure
Etot = relaxed_structure.get_total_energy()
write_cif("final_structure.cif", relaxed_structure)