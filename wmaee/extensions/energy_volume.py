import os
import numpy as np
from wmaee.wmaee import ase_to_pymatgen, working_directory as folder, full_run, VASPInput, parse_output
from wmaee.extensions.elasticity import apply_strain
from wmaee.extensions.common import tqdm
from pymatgen.io.vasp import Incar, Kpoints
from ase import Atoms
from os.path import join, exists




def calc_energy_volume_curve(name, structure, deltas=None, incar=None, kpoints=None, show_output=False, cpus=2,
                             xc_func='gga'):
    if isinstance(structure, Atoms):
        structure = ase_to_pymatgen(structure)

    finished = lambda : exists('vasprun.xml')
    if incar is None:
        incar = Incar(dict(
            NSW=0,
            LREAL='Auto',
            LCHARG='False',
            LWAVE='False'
        ))

    if kpoints is None:
        kpoints = Kpoints.automatic_density_by_vol(32000)

    if deltas is None:
        deltas = np.linspace(-0.1, 0.1, 10)

    hydrostatic_strain = lambda delta: np.eye(3) * delta
    folder_name = lambda seq: 'strain-%i' % seq

    with folder(name):
        data = []
        for i, delta in enumerate(tqdm(deltas, desc='Progress')):
            with folder(folder_name(i)):
                # speed up the calculation by symliniking WAVECAR
                if i == 0:
                    incar['LWAVE'] = True
                else:
                    incar['LWAVE'] = False
                    # make symlink
                    os.system('ln -s %s %s' % (join('..', folder_name(0), 'WAVECAR'), 'WAVECAR'))
                deformed = apply_strain(structure, hydrostatic_strain(delta))
                # print(deformed.lattice.matrix/4.05)
                energy = full_run(VASPInput(incar, deformed, kpoints, xc_func=xc_func), cpus=cpus,
                                  show_output=show_output).final_energy if not finished() else parse_output().final_energy
                data.append((deformed.lattice.volume, energy))
        return list(zip(*data))