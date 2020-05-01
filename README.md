# wmaee
A lightweight tool for the exercises **Werkstoffmodellierung auf atomarer Ebene**.
Provides routines for writing VASP input files, running VASP and parsing its output.
In order to use this module one needs the following modules

* `pymatgen`
* `ase`
* `plotly` 
* `pyiron`


# Installation
Please read through this installation guide carefully. As this module is our own and not available on `pip` or `conda-forge` you have to go through some steps in order to make it run

## General installation
In order to use the convienience function which are shipped with this module all you have to do is to clone the repo. Therefore just execute
```bash
git clone https://git.unileoben.ac.at/cms/wmaee.git
```
in your command line. This will create a new directory with the name **wmaee** in your current folder. Now to use the code you have to add it to your `PYTHONPATH` environment variable

  1. **Permanently**: Add `export PYTHONPATH=/home/dominik/.local/python/wmaee:${PYTHONPATH}` to your `.bashrc` file (This is my example configuration)
  2. Set the environment variable manually whenever you want to use the module

## Make VASP run on your computer with the `wmaee` module
In order to make convenience functions (e.g automatic generation of POTCAR files) work a few steps have to be followed, and a litle bit of configuring is neccessary.

1. Create a configuration directory (e. g `.config` - I have name it so, but you can as you want ) and remember its **absolute** path.
2. Create a directory for each XC funtional and copy VASP's potential archives in there, directly as offered from the vasp webpage
    * In my case I created two folder **gga** and **lda**
3. Have a look in the previously cloned repository there you'll find a file named `defaults.json`
    * **Copy** it into the your **previously** created **configuration directory**
    * `defaults.json` tell `wmaee` which POTCAR file should be taken by default for which species
    * You can simply add more functionals, it's structure is self explanatory
    * **NOTE**: The keys in the `defaults.json` file must match the names of the directories you've chosen to store your POTCAR archives
    * The default keys implemented in the current state are `gga` and `lda`
4. I strongly encourage you to create a folder structure like this shown below

```bash
.config/
├── defaults.json
├── application.json
├── gga
│   └── potpaw_PBE.54.tar.gz
└── lda
    └── potpaw_LDA.54.tar.gz
```
5. Tell the `wmaee` module how to call VASP
    * In the repository directory you'll also find in addition to `defaults.json` a file named `application.json`
    * **Copy** it into the your **previously** created **configuration directory**
    * `application.json` tells the `wmaee` module how to call VASP, however it is more or less a summary of all our *run.sh* scripts
    * The file is structured in the following way, which enables us to resuse it on all our HPC clusters.
    * Now search for the part **"local"**, which represents the configuration of the local machine
    * **"default"** represents the partition, since you probably have no SLURM running on you local machine that is what you want to edit
    * Adapt the configuration so that it fits to your local machine
        * The file on the repo contains the configuration for my workstation -> that gives you a starting point
        * **"preamble"** is a list of commands which are executed before. E. g. load modules on HPC clusters
        * **"binary"** is the binary itself
        * **"command"** the comand how to execute it. You can use also `mpiexec` or `time`. `{cores}` and `{binary}` are wildcards, since the numbers of CPU's can be set in the code
    * Repeat this step for the **vasp_ncl** and **vasp_gam** section if you want to use them

```json
{
    "vasp_std": {
        "mul-hpc": {
            "Phi": {
                ...
            }
        },
        "smmpmech": {
            "default": 
            {

            }
        },
        "local": {
            "default": {
                "preamble": ['source some/intel.sh'],
                "binary": "/path/to/vasp/vasp_std",
                "command": "mpirun -np {cores} {binary}"
            }
        }
    }
}
```
    
6. Set environment variables. You have to tell `wmaee` where it finds the configuration.
    * Either add `export WMAEE_CONFIG_DIR=/path/to/your/.config` into `.bashrc` or add it manually before every usage
7. By default if not specified otherwise `wmaee` will use `local` configuration section in `application.json` as default, and there it will use the `default` configuration by default, however if you want to be explicit which is always a good practice also set manually or in `.bashrc`
    * `export WMAEE_HOSTNAME=local`
    * `export WMAEE_PARTITION=default`

8. **Summary**: After setting up you configuration directory and modifiying `application.json`, it might be a good idea to append the following four lines to your `.bashrc`

```bash
export PYTHONPATH=/path/to/your/wmaee/repo:${PYTHONPATH}
export WMAEE_CONFIG_DIR=/path/to/your/.config
export WMAEE_HOSTNAME=local
export WMAEE_PARTITION=default
```

On my workstation it looks like this
```bash
export PYTHONPATH=/home/dominik/.local/python/wmaee:${PYTHONPATH}
export WMAEE_CONFIG_DIR=/home/dominik/.local/python/wmaee/.config
export WMAEE_HOSTNAME=local
export WMAEE_PARTITION=default
```

## Make LAMMPS run on your computer with the `wmaee` module
**NOTE:** This will only work in an **Anaconda** environment, since the LAMMPS part is built around `pyiron` which is available on `conda-forge` only (or at least the most important packages it needs to run nicely)

### prerequisites
As said `pyiron` is required to run LAMMPS a few packages are needed. You can run
```bash
conda install -c conda-forge pyiron lammps ovito nglview sqsgenerator
```

