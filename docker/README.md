
# How to create a docker image for the `wmaee` exercises

This is a short guide on how to set up the image. The most important file is the `Dockerfile`, however, there will be 
a few more needed. Those are

 - `users.json`: defines users and admins as well as their passwords. **Do never check the `users.json` into the git repo!!! Always use the crypted version** 
 - `crypt.sh`: helper utility to crypt the users JSON file.
 - `uncrypt.sh`: helper utility to uncrypt files.
 - `passwd`: a file containing the cryting password in clear text. (`passwd` is required for `crypt` and `uncrypt`)
 - `env.yaml`: the specifications for the conda env, in which the students will run their calculations
 - `configure_users.py`: a simple python helper script that executes system commands, creates the users and sets up directories, for the users
 - `jupyterhub_config.py`: the configuration file for the JupyterHub
 - `jupyterhub_config_rtc.py`: the real-time collaboration (RTC) configuration file for the JupyterHub (is concatenated with `jupyterhub_config.py`)
 - `wmaee.conf.yaml`: the configuration file for the *wmaee* module

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
located in the same directory as the `Dockerfile` when building the image. In case you change the filenames  you have to adapt the `Dockerfile` accordingly. For safety reason there are copies of the archives on the *fileserver* 
(*dgehringer@fileserver:/share/homes/dgehringer/calculations/wmaee*).
In case you change the paths within the PAW potentials do not forget to adapt the paths in `wmaee.conf.yaml` as well.



## Guide

Please follow the steps closely, and in the correct order! 
**If not specified otherwise all instruction are to be carried out in `wmaee/docker`**

### Prepare the configuration files
Many of the configuration files contain sensitve data, e.g. SSH keys. It is not a good practice to leave the uncrypted on your disk. Therefore, the `Dockerfile` uses crypted versions of these files. 

To crypt files, there is the `crypt.sh` utility. This script reads a password (in clear text) from a file `passwd` and crypts, a file. `uncrypt.sh` does the opposite.

```bash
# will read the passwd file, crypt users.json and output it to users.json.gpg
./crypt.sh users.json
# will read the passwd file, uncrypt users.json.gpg and output it to users.json
./uncrypt.sh users.json.gpg
```

So before starting to build the image creae a file `passwd` which contans a file 

#### Setup the users

The users are configured with a JSON file, named *users.json*. It contains the usernames and the corresponding passwords. Admin users are allowed to shutdown and start servers of other users too. Furthermore, they can log into servers of other non admin users and use to use the RTC feature .

```json
{
  "admins": {
    "dgehringer": "dgehringer",
    "dholec":  "dholec"
  },
  "users": {
    "tleiner": "tleiner",
    "asakic": "asakic"
  }
}
```

This the user file has to be changed each year.

Once all the students and supporting PhD have an entry crypt the file 
```bash
./crypt.sh users.json
```

#### Setup up the access to the fileserver
At the startup of the JupyterHub container it mounts a backup drive via SSHFS to `/media/backup`. To mount a drive via SSHFS, the ssh config must be provided.  For know we back it up onto our fileserver. Thats will now work with any server and user. Four ingredients are needed.

1. a private key
	1. copy your private key into the current directory (*wmaee/docker*) and name it `id_rsa`
	2. crypt it with `./crypt.sh id_rsa`
2. the public key of the server (to prevent interactive questions)
	1. extract the fingerpring of the backup server (in our case it is the fileserver) and store it in a file named `known_hosts`
      1. To keep things short and sweet you can copy the *known_hosts* file from your *.ssh* if you have an access to that server
      2. crypt it with `./crypt.sh known_hosts`
3. the location on the backup server
      1. adapt the environment variable "*BACKUP_REMOTE_DIRECTORY*" in the `Dockerfile`, to a location on the backup server where you have access to (e.g `ENV BACKUP_REMOTE_DIRECTORY=/share/homes/dgehringer/calculations/wmaee/data/2023`)

4. a SSH config file
	1. create a config file that defines the remote
	2. **Note:** The `Dockerfile` renames your `id_rsa` file into `id_rsa_fileserver`
	3. A configuration might look like this 
	```ssh
	Host fileserver
	    HostName <the-ip-of-our-fileserver>
	    User <your-username-on-the-fileserver>
	    Port <the-ssh-port-of-our-fileserver>
	    IdentityFile /root/.ssh/id_rsa_fileserver
	```
	4. Note that the `Host`property must match the line `ENV BACKUP_REMOTE=fileserver` in the `Dockerfile`
	5. crypt it with `./crypt.sh config`


#### Obtain the pseudo potentials
In a last step you have to copy the archive with the pseudo potential named `abinit-potentials.tar.gz` and `gpaw-setups-0.9.20000.tar.gz` into *wmaee/docker*


### building the JupyterHub image
build the container, and hope that everything works out, that might take some time

```bash
#execute this in the directory where the Dockerfile is located
sudo docker build -t wmaee .
```

### prepare the configuration files for the *nginx* reverse proxy
Again the nginx proxy runs in another Docker container
**If not specified otherwise all instruction are to be carried out in `wmaee/docker/nginx-reverse`**

1. obtain the SSL key
	1. the `Dockerfile` assumes that it is named `modelling.unileoben.ac.at.key`: copy it into the directory (`wmaee/docker/nginx-reverse`)
2. obtain the SSL certificate
	1. the `Dockerfile` assumes that it is named `modelling.unileoben.ac.at.cert.pem`: copy it into the directory
3. a file with a large prime number
	1. generate it with ` openssl dhparam -out dhparam.pem 2048` 

### building the JupyterHub image
build the container, and hope that everything works out, that might take some time

```bash
#execute this in the directory where the Dockerfile of the reverse proxy image is located
sudo docker build -t nginx-reverse .
```

### run the setup

```bash
# this simple script is located in wmaee/docker
sudo ./start.sh
```