#!bin/bash

cd wmaee

modules="
codes/vasp
codes/vasp/vasp_ml.py
codes/lammps
core/io.py
core/config.py
core/data_structs.py
core/utils.py
core/visualize.py
scopes/eos.py
scopes/md.py
scopes/dos.py
scopes/structural_analysis.py 
scopes/dos.py
utils/
"

pdoc $modules --docformat numpy  -o ../docs
