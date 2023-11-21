from wmaee.core.io import working_directory
from wmaee.codes.lammps_potential import *
from wmaee.codes.lammps_runner import run_lammps
import pandas as pd
import os

print(get_models())
pot_root, pots = get_potentials('eam/alloy', species=['Al'])

pots = pd.DataFrame(pots).sort_values(by='elements', ignore_index='True')
print(pots)

pot = os.path.join(pot_root, pots['pot_file'][1])
print(pot)

lmp_in = f"""
# LAMMPS input script for FCC Al 3x3x3 supercell NVT simulation

# Initialization
units metal
dimension 3
boundary p p p
atom_style atomic

# Define the lattice
lattice fcc 4.05
region box block 0 6 0 6 0 6
create_box 1 box
create_atoms 1 box

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
run 10000
"""

run_lammps(directory='lmp_test', lmp_in=lmp_in)