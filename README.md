# wmaee

A new generation of `wmaee`, currently at version 2.3!

`wmaee` is a collection of modelling and evaluation tools developed by the MUL-CMS group in Leoben.

Online documentation of functions and modules: [https://wmaee.readthedocs.io/en/latest/](https://wmaee.readthedocs.io/en/latest/)

## Installation

1. Basic package, which relies on `numpy`, `ase`, `matplotlib`, and `scipy`:
   ```bash
   pip install git+https://github.com/MUL-CMS/wmaee.git@master
   ```
   
   You may want/need to force an uninstall first:
   ```bash
   pip uninstall wmaee
   ```

2. Optional functionality via pip extras:
    - atomistic simulations (`lammps`, `kimpy`, `kim-query`):
      ```bash
      pip install wmaee[atomistics]
      ```
    - DFT modules (`gpaw`):
      ```bash
      pip install wmaee[dft]
      ```
   - additional functionality for pyiron (`pyiron_atomistics`, `pyiron_base`):
     ```bash
     pip install wmaee[pyiron]
     ```
    External software such as *VASP* may be required for some modules to work.


## Configuration (`wmaee.conf.yaml`)

TO BE EXTENDED/DOCUMENTED

- `location`: paths to data and resources
- `envpath`: environment paths for executables
- `content`: list of available modules and scripts


## Examples

TO DO
- link to example notebooks

## Contribute

Branch off your own branch, develop and test your own functionality. Once working,request a pull into the main branch.

### Version History

#### v2.3 (2025-09-29)
- `codes/pyiron/pyiron_GRACE_job.py` - implementation of GRACE UMLIPs for pyiron workflows
- update of docs
- slightly polished repo

#### v2.2 (2025-spring)
 - update of $C_{ij}$ projections in `scopes/cij.py`

#### v2.1 (2025-winter)
- `codes/pyiron/pyiron_CHGNet_job.py` - implementation of CHGNet UMLIPs for pyiron workflows
- `codes/pyiron/pyiron_NEB_job.py` - implementation of Nudged Elastic Band (NEB) calculations with VASP for pyiron workflows

#### v2.0 (2024)
- major update based on previous experience, removal of core dependence on pyiron to make it a lighter, more stand alone package.

#### v1.0 (2020-2023)
- initial version initiated by Dominik Gehringer

### 101 on merging branches

To merge branch `david` into `wmaee2`:

```bash
git checkout wmaee2     # switch to branch wmaee2
git merge david         # merge david into wmaee2
git push                # sync local changes
git checkout david      # switch back to david
git rebase wmaee2       # update david with latest wmaee2 changes
```