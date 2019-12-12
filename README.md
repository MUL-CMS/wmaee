# wmaee
A lightweight tool for the exercises **Werkstoffmodellierung auf atomarer Ebene**.
Provides routines for writing VASP input files, running VASP and parsing its output.
In order to use this module one needs the following modules

* `pymatgen`
* `ase`
* `plotly` 

Currently to keep the module lightweight the configuration is hardcoded for the `smmpmech` cluster
## Installation
Plese read through this installation guide carefully
### Obtain the the configuration files
In order to make convenience functions (e.g automatic generation of POTCAR files) work a few steps have to be followed

* Create a configuration directory (e. g `.config`) and remember its **absolute** path.
* Create a directory for each XC funtional and copy VASP's potential archives in there
* Place the file `defaults.json` in the configuration directory. You configuration directory should look something like 
this. The file is available on `smmpmech` cluster (`/calc/dnoeger/software/defaults.json`)
```bash
.config/
├── defaults.json
├── application.json
├── gga
│   └── potpaw_PBE.54.tar.gz
└── lda
    └── potpaw_LDA.54.tar.gz
```
### Adapt the module
After cloning this repository please , modify the function `_get_configuration_driectory` in such a way that it 
returns the **absolute** path of the previously created configuration directory. 
```python
# This is a sample configuration, you have to adapt it to your purposes
def _get_configuration_directory():
    directory_name = '.config'
    prefix = '/calc/dnoeger/cms-exercise'
    return join(directory_name, prefix)
```

### Possibility to adapt it also for `mul-hpc`
In order to use it on the new cluster, also go through the installation procedure from above.
There is already the configuration added for the `mul-hpc` cluster, however it is commented out.
Thus you just have the uncomment the VASP commands in the header of the `wmaae/wmaee.py` file right after the imports.
```python
#__VASP_PREAMBLE = [
#    'module purge',
#    'module load intel',
#    'module load mvapich2/2.2',
#    'module load mkl',
#    'module load scalapack/2.0.2',
#    'ulimit -s unlimited'
#]
#__VASP_COMMAND = 'mpirun -np {tasks} /calc/dnoeger/software/vasp-intel-mvapich2-mkl/5.4.1/bin/vasp_std'
```

### Add the module to your `PYTHONPATH` variable
In order to make the module available to your installation add the path to the cloned repository to your 
`PYTHONPATH` variable either directly in your `.bashrc` or explicitly before you start a Jupyter/IPython or python 
instance.
```bash
export PYTHONPATH="$PYTHONPATH:/path/to/your/cloned/wmaee/repository"
```