from pyiron_atomistics.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistics import pyiron_to_ase, ase_to_pyiron
from pyiron_base import GenericParameters
from ase.io import write, read
from os.path import join
import pandas as pd
from ase import Atoms
from ase.io.trajectory import Trajectory


class CHGNet(AtomisticGenericJob):
    def __init__(self, project, job_name):
        super(CHGNet, self).__init__(project, job_name)
        self.structure = None
        self.input = CalcParams()
        self._executable_activate(codename='chgnet')

    @property
    def fmax(self):
        return self.input.get('fmax')

    @fmax.setter
    def fmax(self, fmax: float):
        self.input.set(fmax=fmax)

    # I use max_iter (because pyiron uses same paramater) to specify max. number of steps during relaxtaion (Default=500)
    # However, in max_iter will be written as steps into .py file, such that chgnet can read it
    # In order to specify type: job.max_iter = 1000
    @property
    def max_iter(self):
        return self.input.get('steps')

    @max_iter.setter
    def max_iter(self, max_iter: int):
        self.input.set(steps=max_iter)

    def write_input(self):
        write(
            join(self.working_directory, 'structure.cif'),
            pyiron_to_ase(self.structure)
        )
        if self._generic_input['calc_mode'] == 'static':
            self._write_calc_static()

    def _write_calc_static(self):
        params = []
        for p in self.input.keys():
            params.append(f'{p}={self.input.get(p)}')
        params = ','.join(params)
        script = [
            'from chgnet.model import CHGNet, StructOptimizer',
            'from ase.io import read',
            'from ase.io.cif import write_cif',
            'from ase.io.trajectory import Trajectory',
            'from pymatgen.io.ase import AseAtomsAdaptor',
            'struct = read("structure.cif")',
            'relaxer = StructOptimizer()',
            'result = relaxer.relax(struct, save_path="relax.traj",'+params+')',
            'final_structure = result["final_structure"]',
            'ase_structure = AseAtomsAdaptor.get_atoms(final_structure)',
            'write_cif("final_structure.cif", ase_structure)'
        ]  # save_path="relax.out" saves the trajectory, how to read and parse the traj file?
        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

    def _parse_calc_static(self, output='error.out'):
        try:
            header = pd.read_csv(
                join(self.working_directory, output),
                skiprows=2,
                nrows=0,
                sep='\s+',
            )
            df = pd.read_csv(
                join(self.working_directory, output),
                skiprows=3,
                sep='\s+',
                header=None
            )
            df.columns = ['optimizer_class'] + list(header.columns[:])
            return df
        except:
            self.logger.warning(f'Cannot partse output from {output} file.')

    def collect_output(self):
        if self._generic_input['calc_mode'] == 'static':
            df = self._parse_calc_static()

            # Parse the relaxed structure (here should be the parsing of the trajectory file)
            # save relaxed structure as cif file
            relaxed_ase_structure = read(
                join(self.working_directory, "final_structure.cif"), format="cif")
            relaxed_pyiron_structure = ase_to_pyiron(relaxed_ase_structure)

            with self.project_hdf5.open('output/generic') as h5out:
                # needs modification, because cells appears as directory in output/generic
                h5out['cells'] = relaxed_pyiron_structure.get_cell()
                h5out['energy_pot'] = df['Energy'].values
                h5out['energy_tot'] = df['Energy'].values
                h5out['max_force'] = df['fmax'].values
                h5out['steps'] = df["Step"].values
                h5out['positions'] = relaxed_pyiron_structure.positions

            # Store the relaxed structure in the HDF5 group
            with self.project_hdf5.open('output') as h5out:
                # Pass the h5out group explicitly
                relaxed_pyiron_structure.to_hdf(h5out)

            self.structure = relaxed_pyiron_structure  # Update the job's structure
            # self.output.cells = relaxed_pyiron_structure.get_cell()
            self.status.finished = True

    def get_structure(self, frame=-1):
        # needs rewrite, because i only read out structure details from the relaxed structure
        num_structures = self.number_of_structures
        if num_structures == 0:
            raise ValueError("No structures found in the job output.")

        # Ensure frame is valid
        if frame < 0:
            frame += num_structures
        if frame >= num_structures:
            raise IndexError(
                f"Frame {frame} is out of range for {num_structures} structures.")

        # Read positions and cell
        with self.project_hdf5.open('output/structure') as h5out:
            positions = h5out["positions"]  # correct
            cell = h5out["cell/cell"]  # correct
            species = h5out["species"]
            indices = h5out["indices"]

        # Generate symbols list
        symbols = [species[i] for i in indices]
        return ase_to_pyiron(Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True))


class CalcParams(GenericParameters):
    def __init__(self, input_file_name=None, table_name="calc_params"):
        super(CalcParams, self).__init__(
            table_name=table_name,
        )
