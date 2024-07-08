# import numpy as np
from wmaee.codes.vasp.vasp_ml import generate_ML_AB
from wmaee.core.io import working_directory

with working_directory('/home/david/work/2024_02_Andritz/data/raw/C48Fe144_slab'):
    generate_ML_AB(overwrite=True)