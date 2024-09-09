# wmaee

A new generation of `wmaee`, version 2.1

Online documentation of functions and modules: [https://wmaee.readthedocs.io/en/latest/](https://wmaee.readthedocs.io/en/latest/)

## wmaee.conf.yaml
- location
- envpath
- content

## Getting the Python environment

1. Install [miniforge](https://github.com/conda-forge/miniforge#install) with mamba solver
2. Install recommended packages: `mamba env update --file env.yaml`  <br />
(This will not create a new conda environment, but instead update the current one! To install it into a new environment, use `mamba env --name wmaee --file env.yaml`)

## Basics of git

### Merging

To merge branch `david` into `wmaee2`:

1. `git checkout wmaee2`: switch to branch `wmaee2`
2. `git merge david`: merge david into `wmaee2`
3. `git push`: sync local changes with the server
4. `git checkout david`: switch back to `david` for further developments
5. `git rebase wmaee2`: sync branch `david` with the latest `wmaee2` content
