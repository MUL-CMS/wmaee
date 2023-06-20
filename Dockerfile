
FROM jupyter/base-notebook

LABEL maintainer="Dominik Gehringer <dominik.gehringer@unileoben.ac.at>"

USER root

WORKDIR /tmp

RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
    git \
    xclip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN mamba install --yes -c conda-forge \
    'abinit' \
    'numpy' \
    'scipy' \
    'matplotlib' \
    'ase' \
    'gpaw' \
    'lammps' \
    'openkim-models' \
    'frozendict' \
    'kim-query' \
    'kimpy' \
    'nglview' \
    'tqdm' \
    'sympy' && \
    mamba clean --all -f -y


RUN pip install git+https://git.unileoben.ac.at/cms/wmaee.git@code-interfaces
