from tensorpotential.calculator import grace_fm
from ase.io import read
from ase.io.cif import write_cif

initial_structure = read("structure.cif")
initial_structure.calc = grace_fm("GRACE-2L-MP-r5")

# perform the static calculation
energy_pot = initial_structure.get_total_energy()
forces = initial_structure.get_forces()
stresses = initial_structure.get_stress()

# write the results to an output file
with open("log.out", "w") as f:
    f.write("step	energy_pot	forces	stresses\n")
    f.write(f"1\t{energy_pot}\t{forces.tolist()}\t{stresses.tolist()}\n")

 # write the final structure to a CIF file
write_cif("final_structure.cif", initial_structure)