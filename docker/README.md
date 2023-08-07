
# How to create a docker image for the `wmaee` exercises

This is a short guide on how to set up the image. The most important file is the `Dockerfile`, however, there will be 
a few more needed. Those are

 - `users.json`: defines users and admins as well as their passwords. **Do never check the `users.json` into the git repo!!! Always use the crypted version** 
 - `crypt.sh`: helper utility to crypt the users JSON file.
 - `uncrypt.sh`: helper utility to uncrypt files.
 - `passwd`: a file containing the cryting password in clear text. (`passwd` is required for `crypt` and `uncrypt`)
 - `env.yaml`: the specifications for the conda env, in which the students will run their calculations
 - `configure_users.py`: a simple python helper script that executes system commands, creates the users and sets up directories
 - `jupyterhub_config.py`: the configuration file for the jupyter hub
 - `wmaee.conf.yaml`: the config file for the *wmaee* module

## Prerequisites

### System requirements

The whole deployment is based on [Docker](https://en.wikipedia.org/wiki/Docker_(software)) while our specific deployment
is based on the [Jupyterhub](https://hub.docker.com/r/jupyterhub/jupyterhub/) image. To create the single user
JupyterLab instances, we use the [systemd](https://github.com/jupyterhub/systemdspawner) spawner. The reasons therefore are:

  - limit CPUs per user
  - limit memory per user
  - isolate user directories
  - hide system devices
  - that leads to improved security

However, this choice comes at a cost. While the userland comes from the image the Docker engine uses the host systems 
kernel. To make systemd spawner work properly *systemd >= 245* is required. You can check the version with 
`systemctl --verison`. [Here](https://github.com/jupyterhub/systemdspawner#systemd-and-linux-distributions) is a list of 
suggested (host) systems.

### wmaee code itself
The current version of the wmaee module that reuses ASE's interfaces is located in the *code-interfaces* branch.
In case you merge the changes into a different branch you have to adapt the `Dockerfile`. Upon building the image
the latest version is installed via pip and git.

### PAW potentials
to make the codes work smoothly for the students, we are going to need PAW potentials for both ABINIT and GPAW.
Currently, the ABINIT potentials are located in an archive named `abinit-potentials.tar.gz`, but you can download them 
from the [ABINIT site of ase](https://wiki.fysik.dtu.dk/ase/ase/calculators/abinit.html). Similarly, the 
GPAW potentials `gpaw-setups-0.9.20000.tar.gz` might be downloaded from the GPAW website. Make sure that those are 
located in the same directory as the `Dockerfile` when building the image. In case the filenames change you have to adapt
the `Dockerfile` accordingly. For safety reason there are copies of the archives on the fileserver 
(*dgehringer@fileserver:/share/homes/dgehringer/calculations/wmaee*).
In case you change the paths within the PAW potentials do not forget to adapt the paths in `wmaee.conf.yaml` as well.



## Guide

### building the container
To setup the container, the scripts need a JSON file specifying the users. More precisely an encrypted version of it.
Such a file might look like: 

```json
{
  "admins": {
    "dgehringer": "1625acf3b1f24be9b36d6ad95ac6c15e",
    "dholec":  "5b066a6f119042feb460f948f73d70c1"
  },
  "users": {
    "tleiner": "08ca0637ae4e440c803367f0a915a6fe"
  }
}
```

once you have gathered all files, run `docker build -t wmaee .` in this directory. This may take a while