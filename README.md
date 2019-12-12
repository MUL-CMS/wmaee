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
### Add the module to your `PYTHONPATH` variable
In order to make the module available to your installation add the path to the cloned repository to your 
`PYTHONPATH` variable either directly in your `.bashrc` or explicitly before you start a Jupyter/IPython or python 
instance.
```bash
export PYTHONPATH="$PYTHONPATH:/path/to/your/cloned/wmaee/repository"
```
