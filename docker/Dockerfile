
FROM jupyterhub/jupyterhub

LABEL maintainer="Dominik Gehringer <dominik.gehringer@unileoben.ac.at>"

ENV REPOSITORY_URL=https://git.unileoben.ac.at/cms/wmaee.git
ENV REPOSITORY_BRANCH=code-interfaces

ENV JUPYTERHUB_CONF_DIR=/etc/jupyterhub
ENV JUPYTERHUB_BASE_DIR=/opt/jupyterhub
ENV JUPYTERHUB_USERS_DIR=/opt/jupyterhub/users

ENV MAMBA_URL=https://micro.mamba.pm/api/micromamba/linux-64/latest
ENV MAMBA_BIN=/bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/mamba

USER root

RUN mkdir -p $JUPYTERHUB_BASE_DIR
RUN mkdir -p $JUPYTERHUB_USERS_DIR
RUN mkdir -p $JUPYTERHUB_CONF_DIR

WORKDIR /

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    git \
    xclip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install micromamba
RUN mkdir -p $MAMBA_ROOT_PREFIX
RUN bash -c "$(curl https://micromamba.pfx.dev/install.sh)"

RUN mv ~/.local/bin/micromamba /bin/micromamba

RUN micromamba config append channels conda-forge
RUN micromamba self-update

WORKDIR $JUPYTERHUB_CONF_DIR

#install systemd spawner
RUN pip install jupyterhub-systemdspawner
RUN pip install jupyterlab

# create jupyter hub
COPY env.yaml .

RUN micromamba create -f env.yaml
RUN bash -c "eval $(micromamba shell hook --shell bash) && micromamba activate wmaee && pip install git+$REPOSITORY_URL@$REPOSITORY_BRANCH"

# publish the mamba environment, by installing it as a foreign kernel
#https://ipython.readthedocs.io/en/latest/install/kernel_install.html#kernels-for-different-environments
RUN bash -c "eval $(micromamba shell hook --shell bash) && micromamba activate wmaee && python -m ipykernel install --prefix=$MAMBA_ROOT_PREFIX/envs/wmaee --name wmaee --display-name \"Python 3.8 (wmaee)\""

ENV POTENTIAL_DIR_ABINIT=/opt/potentials/abinit
ENV POTENTIAL_DIR_GPAW=/opt/potentials/gpaw

RUN mkdir -p $POTENTIAL_DIR_ABINIT
RUN mkdir -p $POTENTIAL_DIR_GPAW

WORKDIR $POTENTIAL_DIR_ABINIT
COPY abinit-potentials.tar.bz2 .
RUN tar xf abinit-potentials.tar.bz2


WORKDIR $JUPYTERHUB_CONF_DIR
# create system users
COPY users.json.gpg .
COPY uncrypt.sh .
COPY passwd .
COPY configure_users.py .
COPY jupyterhub_config.py .

# encrypt the file
RUN chmod +x uncrypt.sh
RUN bash -c "./uncrypt.sh users.json.gpg"

# set up the users
RUN python configure_users.py users.json jupyterhub_config.py
RUN rm passwd configure_users.py users.json users.json.gpg uncrypt.sh




CMD ["jupyterhub", "-f", "/etc/jupyterhub/jupyterhub_config.py"]


