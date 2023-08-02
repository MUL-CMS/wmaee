# we can pass function to the dictionary as values
# https://jupyterhub-tutorial.readthedocs.io/en/latest/spawners.html
import os

MAMBA_ROOT_PREFIX="/opt/mamba"
ENV_NAME = "wmaee"

c.Spawner.environment = {
    "PATH": os.path.join(MAMBA_ROOT_PREFIX, 'condabin') + ":" + os.path.join(MAMBA_ROOT_PREFIX, "envs", ENV_NAME, "bin") + ":" + os.environ.get("PATH"),
    "MAMBA_ROOT_PREFIX": MAMBA_ROOT_PREFIX,
    "WMAEE_CONFIG_FILE": os.path.join("/etc", ENV_NAME, "wmaee.conf.yaml"),
    "JUPYTER_DATA_DIR": lambda spawner: f"/opt/jupyterhub/users/{spawner.user.name}/.local/share/jupyter",
    "JUPYTER_RUNTIME_DIR": lambda spawner: f"/opt/jupyterhub/users/{spawner.user.name}/.local/share/jupyter"
}

c.JupyterHub.authenticator_class = 'jupyterhub.auth.PAMAuthenticator'
# systemd spawner configuration
c.SystemdSpawner.dynamic_users = False
c.SystemdSpawner.mem_limit = '5G'
c.SystemdSpawner.cpu_limit = 2
c.SystemdSpawner.user_workingdir = '/opt/jupyterhub/users/{USERNAME}'
c.SystemdSpawner.username_template = '{USERNAME}'
c.SystemdSpawner.default_shell = '/bin/bash'
c.SystemdSpawner.isolate_devices = True
c.SystemdSpawner.disable_user_sudo = True
c.SystemdSpawner.extra_paths = ['/opt/mamba/envs/wmaee/bin']
c.SystemdSpawner.readwrite_paths = ['/opt/jupyterhub/users/{USERNAME}']
c.SystemdSpawner.readonly_paths = ['/']
