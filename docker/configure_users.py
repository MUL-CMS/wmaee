

import os
import sys
import json
from typing import *

JUPYTERHUB_USERS_DIR = os.environ["JUPYTERHUB_USERS_DIR"]
MAMBA_ROOT_PREFIX = os.environ["MAMBA_ROOT_PREFIX"]
BACKUP_DIRECTORY = os.environ["BACKUP_DIRECTORY"]
ENV_NAME = "wmaee"

_, user_file, config_file, *_ = sys.argv


def system(cmd: str, exit_on_fail: bool = True) -> int:
    print(f"[user_config]: \"{cmd}\"")
    result = os.system(cmd)
    if exit_on_fail:
        assert result == 0
    return result


def create_users(u: Dict[str, str]) -> NoReturn:
    for username, password in u.items():
        user_directory = os.path.join(JUPYTERHUB_USERS_DIR, username)
        system(f"useradd -m -d {user_directory} -p $(echo \"{password}\" | openssl passwd -1 -stdin) -s /bin/bash {username}")
        # create kernel directory
        jupyter_directory = os.path.join(user_directory, ".local", "share", "jupyter")
        system(f"mkdir -p {jupyter_directory}")
        # create symlinks to the folders in the main environment
        folders = {"kernels", "lab", "labextensions", "nbconvert", "nbextensions"}
        env_jupyter_base_path = os.path.join(MAMBA_ROOT_PREFIX, "envs", ENV_NAME, "share", "jupyter")
        # create the symlinks to the jupyterlab installation in /opt/mamba env, since we use this with a "foreign" kernelspec
        for folder in folders:
            system(f"cd {jupyter_directory} && ln -s {os.path.join(env_jupyter_base_path, folder)} {folder}")

        # init micromamba shell-hook for each user separately
        system(f"su {username} -c \"/bin/micromamba shell init --shell bash --root-prefix={MAMBA_ROOT_PREFIX}\"")
        # let the user own the .local directory in his home folder
        system(f"chown -R {username}:{username} {os.path.join(user_directory, '.local')}")



def format_user_set(users: Iterable[str]) -> str:
    contents = ", ".join(f"'{user}'" for user in users)
    return f"{{{contents}}}"


if __name__ == "__main__":
    with open(user_file) as handle:
        user_conf = json.load(handle)

    admins = user_conf.get("admins")
    users = user_conf.get("users")

    all_users = admins.copy()
    all_users.update(users)
    create_users(all_users)
    # append the users to the notebook auth
    system(f"echo \"c.Authenticator.admin_users = {format_user_set(admins)}\" >> {config_file}")
    system(f"echo \"c.Authenticator.allowed_users = {format_user_set(all_users)}\" >> {config_file}")





