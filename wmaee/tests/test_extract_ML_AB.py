# import numpy as np
from wmaee.codes.vasp.vasp_ml import generate_ML_AB
from wmaee.core.io import working_directory

with working_directory('/home/david/work/2024_02_Andritz/data/raw/test'):
    generate_ML_AB(input='OUTCAR', input_type='OUTCAR', overwrite=True)
    # generate_ML_AB(input='vasprun.xml', input_type='vasprun', overwrite=True)