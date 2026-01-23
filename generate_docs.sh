#!bin/bash

cd wmaee

modules="
codes/vasp
codes/vasp/vasp_ml.py
codes/lammps
codes/pyiron/pyiron_cluster.py
codes/pyiron/pyiron_NEB_task.py
codes/pyiron/pyiron_CHGNet_job.py
codes/pyiron/pyiron_GRACE_job.py
core/io.py
core/config.py
core/data_structs.py
core/utils.py
core/visualize.py
scopes/cij.py
scopes/debye.py
scopes/eos.py
scopes/md.py
scopes/dos.py
scopes/structural_analysis.py 
scopes/dos.py
utils/
"

pdoc $modules --docformat numpy  -o ../docs
