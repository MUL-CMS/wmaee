from wmaee.core.io import working_directory
import wmaee.codes.lammps as lammps
import pandas as pd
import numpy as np
import os

print(lammps.get_models())
pot_root, pots = lammps.get_potentials('eam/alloy', species=['Al'])

pots = pd.DataFrame(pots).sort_values(by='elements', ignore_index='True')
print(pots)

pot = os.path.join(pot_root, pots['pot_file'][1])
print(pot)

from ase import Atoms
struct = Atoms(
    '4Al',
    cell = 4.05*np.eye(3),
    scaled_positions = [
        [0, 0, 0],
        [0.5, 0.5, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
        ],
)*[6, 6, 6]

data_in = 'struct_in.data'
lmp_in = f"""
# LAMMPS input script for FCC Al 3x3x3 supercell NVT simulation

# Initialization
units metal
dimension 3
boundary p p p
atom_style atomic

# # Define the lattice
# lattice fcc 4.05
# region box block 0 6 0 6 0 6
# create_box 1 box
# create_atoms 1 box
# # OR
# Read data from data file
read_data {data_in} 

# Set mass for aluminum
mass 1 27.0

# Define potential (replace 'your_potential_file' with the actual potential file)
pair_style eam/alloy
pair_coeff * * {pot} Al

# Set temperature and timestep
timestep 0.002 # default for metal units: 0.001
velocity all create 300.0 12345

# Set thermostat for NVT ensemble with automated setting of damping 
# parameter for the Nose-Hoover thermostat
fix 1 all nvt temp 300.0 300.0 $(100.0*dt)

# Output values every 100 timesteps
thermo 100

# Run simulation
run 1000

thermo_style custom step temp pe etotal press vol
run 1000
"""
with working_directory('lmp_test'):
    lammps.write_lammps_data(data_in, struct)
    lammps.run_lammps(lmp_in=lmp_in)
    output = lammps.parse_logfile('log.lammps')
for partial in output:
    print(partial)