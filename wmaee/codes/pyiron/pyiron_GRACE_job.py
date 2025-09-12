# Metadata
__author__ = "David Holec, Amin Sakic"
__maintainer__ = "David Holec"
__email__ = "david.holec@unileoben.ac.at"
__date__ = "2025-08-03"
__version__ = "0.11"
__license__ = "MIT"

from pyiron_atomistics.atomistics.job.atomistic import AtomisticGenericJob
from pyiron_atomistics import pyiron_to_ase, ase_to_pyiron
from pyiron_base import GenericParameters

from tensorpotential.calculator import grace_fm
from tensorpotential.calculator import TPCalculator


from ase.io import write, read
import io
from os.path import join
from ase import Atoms
from ase.constraints import ExpCellFilter
from ase.io.trajectory import Trajectory


import pandas as pd
import numpy as np
from typing import Literal, Optional


class Grace(AtomisticGenericJob):
    """
    Custom job class for using Grace within the pyiron framework.

    This class provides methods to perform static calculations, structure relaxation,
    and molecular dynamics simulations using Grace umlips.

    Attributes
    ----------
    structure : pyiron.atomistics.structure.atoms.Atoms
        Structure associated with the job.
    input : CalcParams
        Input parameters for the Grace calculations.
    """

    def __init__(self, project, job_name):
        """
        Initialize the Grace job.

        Parameters
        ----------
        project : pyiron.Project
            The pyiron project in which the job resides.
        job_name : str
            Name of the job.
        """
        super(Grace, self).__init__(project, job_name)
        self.structure = None
        self.input = CalcParams()
        if self.project_hdf5.file_exists:
            self.input.from_hdf(self.project_hdf5, group_name='input')
        # executable is in pyiron_resources/grace/bin
        self._executable_activate(codename='grace')

    @property
    def fmax(self) -> float:
        """
        Maximum force tolerance for ionic relaxation.

        This property sets the maximum force tolerance (in eV/Å) used during ionic relaxation.

        Returns
        -------
        float
            Maximum force tolerance (eV/Å).
        """
        return self._generic_input['ionic_force_tolerance']

    @fmax.setter
    def fmax(self, fmax: float):
        """
        Set the maximum force tolerance for ionic relaxation.

        Parameters
        ----------
        fmax : float
            Maximum force tolerance (eV/Å).
        """
        self._generic_input['ionic_force_tolerance'] = fmax

    @property
    def max_iter(self) -> int:
        """
        Get the maximum number of iterations for relaxation.

        Returns
        -------
        int
            Maximum number of iterations.
        """
        return self._generic_input['n_print']

    @max_iter.setter
    def max_iter(self, max_iter: int):
        """
        Set the maximum number of iterations for relaxation.

        Parameters
        ----------
        max_iter : int
            Maximum number of iterations.
        """
        self._generic_input['n_print'] = max_iter

    def calc_static(
        self,
        grace_model: Literal["GRACE-1L-MP-r6", "GRACE-2L-MP-r5",
                             "GRACE-2L-MP-r6", "GRACE-FS-OAM",
                             "GRACE-1L-OAM", "GRACE-2L-OAM",
                             "GRACE-FS-OMAT", "GRACE-1L-OMAT",
                             "GRACE-2L-OMAT"] = "GRACE-2L-OMAT",
    ):
        """
        Perform a static calculation.

        Parameters
        ----------
        grace_model : {"GRACE-1L-MP-r6", "GRACE-2L-MP-r5", "GRACE-2L-MP-r6", "GRACE-FS-OAM",
        "GRACE-1L-OAM", "GRACE-2L-OAM", "GRACE-FS-OMAT", "GRACE-1L-OMAT", "GRACE-2L-OMAT"}, optional
            The Grace model to use for the calculation. Default is "GRACE-2L-OMAT".

        This method sets up the Grace calculator with the specified model and prepares
        the job for a static calculation. The structure must be set before calling this method.

        Raises
        ------
        ValueError
            If the structure is not set before calling this method.
        """
        if self.structure is None:
            raise ValueError(
                "Structure must be set before calling calc_static().")

        # Set up the Grace calculator
        calc = grace_fm(grace_model)
        self.structure.calc = calc
        self._generic_input["calc_mode"] = "static"
        self._generic_input["model"] = grace_model

    def _write_calc_static(self):
        """
        Write the calculation script for performing static calculations.

        This method generates a Python script that uses Grace for static calculations
        and writes it to a file in the working directory.
        """
        model = self._generic_input["model"]

        script = [
            'from tensorpotential.calculator import grace_fm',
            'from ase.io import read',
            'from ase.io.cif import write_cif',
            '\ninitial_structure = read("structure.cif")']
        script += [
            f'initial_structure.calc = grace_fm("{model}")',
            '\n# perform the static calculation',
            'energy_pot = initial_structure.get_total_energy()',
            'forces = initial_structure.get_forces()',
            'stresses = initial_structure.get_stress()',
            '\n# write the results to an output file',
            'with open("log.out", "w") as f:',
            '    f.write("step\tenergy_pot\tforces\tstresses\\n")',
            '    f.write(f"1\\t{energy_pot}\\t{forces.tolist()}\\t{stresses.tolist()}\\n")',
            '\n # write the final structure to a CIF file',
            'write_cif("final_structure.cif", initial_structure)',
        ]
        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

    def calc_minimize(
        self,
        grace_model: Literal["GRACE-1L-MP-r6", "GRACE-2L-MP-r5",
                             "GRACE-2L-MP-r6", "GRACE-FS-OAM",
                             "GRACE-1L-OAM", "GRACE-2L-OAM",
                             "GRACE-FS-OMAT", "GRACE-1L-OMAT",
                             "GRACE-2L-OMAT"] = "GRACE-2L-OMAT",
        algorithm: Literal["FIRE", "BFGS", "LBFGS"] = "FIRE",
        algorithm_kwargs: dict = None,
        relax_cell: bool | None = False,
        relax_cell_kwargs: dict = None,
        ionic_force_tolerance: float = 0.1,
        max_iter: int = 500,
        n_print: int = 1,
        save_path: str | None = 'relax.traj'
    ):
        """
        Perform a structure relaxation (geometry optimization).

        Parameters
        ----------
        grace_model : str,
            The Grace model to use for the calculation. Default is "GRACE-2L-OMAT".
        algorithm : str, optional
            Choose an ASE algorithm for the ionic relaxation. Default is "FIRE".
        algorithm_kwargs: dict, optional
            Additional keyword arguments for the chosen algorithm. Default is None.
        relax_cell : bool or None, optional
            Whether to allow relaxation of the simulation cell (default is False).
        relax_cell_kwargs: dict, optional
            Additional keyword arguments for cell relaxation (default is None).
        ionic_force_tolerance : float, optional
            Convergence criterion for forces (default is 0.1 eV/Å).
        max_iter : int, optional
            Maximum number of iterations (default is 500).
        n_print : int, optional
            Frequency of logging output (default is 1).
        save_path : str or None, optional
            Path to save the trajectory file (default is 'relax.traj').
        """
        super().calc_minimize(
            ionic_energy_tolerance=0,
            ionic_force_tolerance=ionic_force_tolerance,
            max_iter=max_iter,
        )

        # Saving the parameters in dictionary
        self._generic_input["model"] = grace_model
        self.input['algo'] = algorithm
        self.input['algo_kwargs'] = algorithm_kwargs or {}
        self.input['relax_cell'] = relax_cell
        self.input['relax_cell_kwargs'] = relax_cell_kwargs or {}
        self._generic_input['ionic_force_tolerance'] = ionic_force_tolerance
        self._generic_input['n_print'] = n_print
        self.input['save_path'] = save_path

    def _write_calc_minimize(self):
        """
        Write the calculation script for performing structure relaxation.

        This method generates a Python script that uses Grace for structure
        relaxation and writes it to a file in the working directory.
        """
        # Extract parameters from dictionary

        # Grace model
        model = self._generic_input["model"]

        # ASE optimizer
        algo = self.input.get('algo', 'FIRE')

        # Print warning if algo not in allowed list
        if algo not in ['FIRE', 'BFGS', 'LBFGS']:
            raise ValueError(
                f"Unsupported optimizer: {algo}. Supported optimizers are: FIRE, BFGS, LBFGS.")

        # ASE optimizer kwargs
        algo_kwargs = self.input.get('algo_kwargs', {}).copy()

        # Add trajectory parameter
        save_path = self.input.get("save_path", "relax.traj")
        print("Trajectory name:", save_path)
        algo_kwargs['trajectory'] = f"{save_path}"

        # Reformat dictionary to string
        algo_kwargs_str = ", ".join(
            f'{k}="{v}"' if isinstance(v, str) else f"{k}={v}"
            for k, v in algo_kwargs.items()
        )

        # Cell relaxation
        relax_cell = self.input.get('relax_cell', False)

        # Cell relaxation kwargs ---> see ase.constraints.ExpCellFilter
        relax_cell_kwargs = self.input.get('relax_cell_kwargs', {})

        # Reformat dictionary to string
        # Convert to string for script
        relax_cell_kwargs_str = ", ".join(
            f'{k}="{v}"' if isinstance(v, str) else f"{k}={v}"
            for k, v in relax_cell_kwargs.items()
        )

        # Set parameters for minimization
        params_run = f"fmax={self._generic_input['ionic_force_tolerance']}, steps={self._generic_input['max_iter']}"

        # Save parameters in
        self._generic_input['calc_mode'] = 'minimize'
        self._generic_input['relax_cell'] = self.input.get('relax_cell', False)
        self._generic_input['save_path'] = self.input.get(
            'save_path', 'relax.traj')

        script = []
        import_statements = ['from tensorpotential.calculator import grace_fm',
                             'from ase.io import read',
                             'from ase.io.cif import write_cif',
                             'from ase.io.trajectory import Trajectory',
                             'from ase.optimize import FIRE, BFGS',
                             'from ase.constraints import ExpCellFilter as ECF',]
        script += import_statements

        init_calc = [f'initial_structure = read("structure.cif")',
                     f'calc = grace_fm("{model}")',
                     'initial_structure.calc = calc',]
        script += init_calc

        if not relax_cell:
            calc_script = [
                f'relaxation = {algo.upper()}(atoms=initial_structure, {algo_kwargs_str}).run({params_run})',]
        else:
            calc_script = [
                f'relaxation = {algo.upper()}(ECF(atoms=initial_structure, {relax_cell_kwargs_str}), {algo_kwargs_str}).run({params_run})',]

        script += calc_script

        parse = ['relaxed_structure = initial_structure',
                 'Etot = relaxed_structure.get_total_energy()',
                 'write_cif("final_structure.cif", relaxed_structure)']
        script += parse

        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

    def calc_md(
        self,
        temperature=300,
        pressure=None,
        n_ionic_steps=100,
        time_step=2.0,  # in femto-seconds
        n_print=10,
        temperature_damping_timescale=100.0,  # taut: float | None = None,
        pressure_damping_timescale=None,  # taup: float | None = None,
        trajectory="md_run.traj",
        logfile="md_run.log",
        on_isolated_atoms: Literal["ignore", "warn", "error"] = "warn",
        starting_temperature: int | None = None,
        ensemble: str = "nvt",
    ):
        """
        Not implemented yet.

        Parameters
        ----------
        temperature : float, optional
            Target temperature of the simulation in Kelvin (default is 300 K).
        pressure : float or None, optional
            Target pressure for the simulation in bar (default is None).
        n_ionic_steps : int, optional
            Number of ionic steps to simulate (default is 100).
        time_step : float, optional
            Time step for the simulation in femtoseconds (default is 2.0 fs).
        n_print : int, optional
            Frequency of printing simulation logs (default is 10).
        temperature_damping_timescale : float, optional
            Timescale for temperature damping in femtoseconds (default is 100.0 fs).
        pressure_damping_timescale : float or None, optional
            Timescale for pressure damping in femtoseconds (default is None).
        trajectory : str, optional
            Name of the output trajectory file (default is 'md_run.traj').
        logfile : str, optional
            Name of the output log file (default is 'md_run.log').
        on_isolated_atoms : {'ignore', 'warn', 'error'}, optional
            Behavior for isolated atoms (default is 'warn').
        starting_temperature : int or None, optional
            Initial temperature of the simulation (default is the target temperature).
        ensemble : str, optional
            Ensemble type, e.g., 'nvt' or 'npt' (default is 'nvt').
        """
        super().calc_md(
            temperature=temperature,
            pressure=pressure,
            n_ionic_steps=n_ionic_steps,
            time_step=time_step,
            n_print=n_print,
            temperature_damping_timescale=temperature_damping_timescale,
            pressure_damping_timescale=pressure_damping_timescale,
        )
        self.input['ensemble'] = ensemble
        self.input['trajectory'] = trajectory
        self.input['logfile'] = logfile
        self.input['on_isolated_atoms'] = on_isolated_atoms
        if starting_temperature is None:
            starting_temperature = temperature
        self.input['starting_temperature'] = starting_temperature

    def _write_calc_md(self):
        """
        Write the calculation script for performing molecular dynamics (MD).

        This method generates a Python script that uses CHGNet for running MD
        simulations and writes it to a file in the working directory.
        """
        params = []
        params += [f"temperature={self._generic_input['temperature']},pressure={self._generic_input['pressure']},timestep={self._generic_input['time_step']}"]
        params += [f"loginterval={self._generic_input['n_print']},taut={self._generic_input['temperature_damping_timescale']},taup={self._generic_input['pressure_damping_timescale']}"]
        allowed_params = ['ensemble', 'thermostat', 'starting_temperature', 'bulk_modulus',
                          'crystal_feas_logfile', 'append_trajectory', 'use_device', 'trajectory', 'logfile', 'on_isolated_atoms']
        for p in self.input.keys():
            v = self.input.get(p)
            if p in allowed_params and not v == None:
                if isinstance(v, str):
                    params.append(p+"=\""+v+"\"")
                else:
                    params.append(f'{p}={self.input.get(p)}')
        params = ','.join(params)
        fix_imports = ''
        fixes = []
        if 'fix_imports' in self.input.keys():
            fix_imports = self.input.get('fix_imports')
        if 'fixes' in self.input.keys():
            f = self.input.get('fixes')
            if isinstance(f, str):
                fixes = [self.input.get('fixes')]
            elif isinstance(f, list):
                fixes = f
            else:
                self.logger.warning(
                    'Fixes are neither string or array of strings. Ignoring')

        script = [
            'from chgnet.model import CHGNet, MolecularDynamics',
            'chgnet = CHGNet.load()',
            'from ase.io import read',
            'from ase.io.cif import write_cif',
            fix_imports,
            'struct = read("structure.cif")']
        for f in fixes:
            script.append(f'struct.set_constraint({f})')
        script += [
            'md = MolecularDynamics(atoms=struct,model=chgnet,'+params+')',
            f'md.run({self._generic_input["n_ionic_steps"]})',
            'traj = read(f"md_run.traj", index=":")',
            'write_cif("final_structure.cif", traj[-1])'
        ]
        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

    def write_input(self) -> None:
        """
        Write input files required for the Grace calculation.

        This function writes the structure file (`structure.cif`), stores the input
        dictionary to the HDF5 format, and delegates further writing depending on the
        calculation mode (MD or static/minimize).

        Returns
        -------
        None
        """
        write(join(self.working_directory, 'structure.cif'),
              pyiron_to_ase(self.structure))
        self.input.to_hdf(self.project_hdf5, group_name='input')
        if self._generic_input['calc_mode'] == 'md':
            self._write_calc_md()
        elif self._generic_input['calc_mode'] == 'static':
            self._write_calc_static()
        else:
            self._write_calc_minimize()

    def _parse_calc_static(self, output: str = 'log.out', skip: int = 0) -> Optional[pd.DataFrame]:
        """
        Parse the output of a static calculation.

        This function reads the static calculation output file and returns the
        corresponding DataFrame containing the optimization step and energy information.

        Parameters
        ----------
        output : str, optional
            The name of the output file to parse, by default 'log.out'.
        skip : int, optional
            The number of rows to skip at the start of the file, by default 3.

        Returns
        -------
        pd.DataFrame or None
            A DataFrame containing the optimization results, or None if parsing fails.
        """
        try:
            header = pd.read_csv(
                join(self.working_directory, output),
                skiprows=skip,
                nrows=0,
                sep=r'\t+'
            )
            df = pd.read_table(
                join(self.working_directory, output),
                skiprows=skip+1,
                sep=r'\t+',
                header=None
            )
            # df.columns = ['optimizer_class'] + list(header.columns[:])
            df.columns = list(header.columns[:])

            print(df.columns)
            return df
        except:
            self.logger.warning(f'Cannot parse output from {output} file.')
            return None

    def _parse_calc_minimize(self, output: str = 'error.out') -> Optional[pd.DataFrame]:
        """
        Parse the output of a minimization.

        This function reads the minimization output file and returns the
        corresponding DataFrame containing the optimization step and energy information.

        Parameters
        ----------
        output : str, optional
            The name of the output file to parse, by default 'error.out'.

        Returns
        -------
        pd.DataFrame or None
            A DataFrame containing the optimization results, or None if parsing fails.
        """
        # ASE optimizer
        algo = self.input.get('algo', 'FIRE')

        try:
            # Find the header line dynamically
            with open(join(self.working_directory, output), 'r') as f:
                lines = f.readlines()
            for i, line in enumerate(lines):
                if line.strip().startswith('Step'):
                    header_idx = i
                    break
            else:
                self.logger.warning(f'No table header found in {output}.')
                return None

            table_lines = []
            # Ensure that the line contains the expected algorithm
            for line in lines[header_idx:]:
                if algo in line.strip() or "Step" in line.strip():
                    table_lines.append(line)

            # Read table into DataFrame
            table_str = ''.join(table_lines)
            df = pd.read_csv(io.StringIO(table_str),
                             sep=r'\s+', engine='python')

            return df
        except Exception as e:
            self.logger.warning(f'Cannot parse output from {output} file: {e}')
            return None

    def _parse_calc_md(self, output: str = 'md_run.log') -> Optional[pd.DataFrame]:
        """
        Parse the output of a molecular dynamics (MD) simulation.

        This function reads the MD log file and returns the corresponding DataFrame
        containing time, total energy, potential energy, kinetic energy, and temperature.

        Parameters
        ----------
        output : str, optional
            The name of the MD log file to parse, by default 'md_run.log'.

        Returns
        -------
        pd.DataFrame or None
            A DataFrame containing the MD simulation results, or None if parsing fails.
        """
        try:
            df = pd.read_csv(
                join(self.working_directory, output),
                skiprows=1,
                sep=r'\s+'
            )
            df.columns = ['Time', 'Etot', 'Epot', 'Ekin', 'T']
            return df
        except:
            self.logger.warning(f'Cannot parse output from {output} file.')
            return None

    def collect_output(self) -> None:
        """
        Collect and store the results of the calculation.

        This function parses the final structure, trajectory, and other output properties,
        and stores them in the HDF5 format. The structure is updated, and various properties
        such as energies, forces, and magnetization are saved to the output group.

        Returns
        -------
        None
        """

        def voigt_to_tensor(s):
            # s is a 6-element array from atoms.get_stress()
            return np.array([
                [s[0], s[5], s[4]],
                [s[5], s[1], s[3]],
                [s[4], s[3], s[2]]
            ])

        relaxed_ase_structure = read(
            join(self.working_directory, 'final_structure.cif'), format='cif')
        relaxed_pyiron_structure = ase_to_pyiron(relaxed_ase_structure)

        # Store the relaxed structure in the HDF5 group
        with self.project_hdf5.open('output') as h5out:
            h5out['structure'] = relaxed_pyiron_structure
        self.structure = relaxed_pyiron_structure  # Update the job's structure

        # Parse trajectory from minimize
        if self._generic_input['calc_mode'] == 'minimize':
            # For minimize, we expect a trajectory file with the relaxed structure
            traj_file = self.input['save_path']
            try:
                file = join(self.working_directory, traj_file)
                # use ASE Trajectory to read the file
                traj = Trajectory(file)
                with self.project_hdf5.open('output/generic') as h5out:
                    h5out['cells'] = np.array(
                        [atoms.get_cell() for atoms in traj])
                    h5out['positions'] = np.array(
                        [atoms.positions for atoms in traj])
                    h5out['stresses'] = np.array(
                        [voigt_to_tensor(atoms.get_stress()) for atoms in traj])
                    h5out['forces'] = np.array(
                        [atoms.get_forces() for atoms in traj])
            except:
                self.logger.warning(
                    f'Cannot parse output from {traj_file} file.')
        elif self._generic_input['calc_mode'] == 'md':
            # Parse trajectory from MD
            traj_file = self.input['trajectory']
            try:
                traj = read(join(self.working_directory, traj_file), index=":")
                with self.project_hdf5.open('output/generic') as h5out:
                    h5out['cells'] = np.array([s.get_cell() for s in traj])
                    h5out['positions'] = np.array([s.positions for s in traj])
                    h5out['stresses'] = np.array(
                        [s.get_stress() for s in traj])
                    h5out['forces'] = np.array([s.get_forces() for s in traj])
            except:
                self.logger.warning(
                    f'Cannot parse output from {traj_file} file.')

        # Parse other properties
        if self._generic_input['calc_mode'] == 'static':
            df = self._parse_calc_static()
            with self.project_hdf5.open('output/generic') as h5out:
                h5out['energy_pot'] = df['energy_pot'].values
                h5out['energy_tot'] = df['energy_pot'].values
                h5out['max_force'] = np.array(df['forces'].values).max()
                h5out['steps'] = df['step'].values
        elif self._generic_input['calc_mode'] == 'minimize':
            df = self._parse_calc_minimize()
            with self.project_hdf5.open('output/generic') as h5out:
                h5out['energy_pot'] = df['Energy'].values
                h5out['energy_tot'] = df['Energy'].values
                h5out['forces'] = df['fmax'].values
                h5out['steps'] = df['Step'].values
        elif self._generic_input['calc_mode'] == 'md':
            df = self._parse_calc_md(self.input['logfile'])
            with self.project_hdf5.open('output/generic') as h5out:
                h5out['energy_pot'] = df['Epot'].values
                h5out['energy_tot'] = df['Etot'].values
                h5out['energy_kin'] = df['Ekin'].values
                h5out['temperature'] = df['T'].values
                h5out['time'] = df['Time'].values
                h5out['steps'] = df.index.values

        self.status.finished = True

    def get_structure(self, frame: int = -1) -> Atoms:
        """
        Retrieve the structure at a given frame from the calculation output.

        This function returns the atomic structure (positions, cell, species) at the
        specified frame index. The frame index is adjusted to ensure it is valid.

        Parameters
        ----------
        frame : int, optional
            The frame index to retrieve. If negative, it is counted from the last frame,
            by default -1.

        Returns
        -------
        Atoms
            The ASE Atoms object representing the structure at the given frame.

        Raises
        ------
        ValueError
            If no structures are found in the job output.
        IndexError
            If the frame index is out of range.
        """
        num_structures = self.number_of_structures
        if num_structures == 0:
            raise ValueError("No structures found in the job output.")

        # Ensure frame is valid
        if frame < 0:
            frame += num_structures
        if frame >= num_structures:
            raise IndexError(
                f"Frame {frame} is out of range for {num_structures} structures.")

        with self.project_hdf5.open('output/generic') as h5out:
            positions = h5out["positions"][frame]
            cell = h5out["cells"][frame]
        # Read positions and cell
        with self.project_hdf5.open('output/structure') as h5out:
            species = h5out["species"]
            indices = h5out["indices"]

        # Generate symbols list
        symbols = [species[i] for i in indices]
        return ase_to_pyiron(Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True))


class CalcParams(GenericParameters):
    """
    Class to manage the calculation parameters for Grace calculations.

    This class extends from `GenericParameters` and is used to define and manage
    parameters specific to Grace calculations. It allows for initialization of
    parameters and stores them in a table, with a default table name of "grace".
    The class can be used to load and store various Grace calculation settings and
    other related configurations.

    Parameters
    ----------
    table_name : str, optional
        The name of the table where the parameters are stored. By default, it is set
        to "grace", but it can be customized if necessary.

    Methods
    -------
    __init__(self, table_name="grace")
        Initializes the parameters for the calculation.
    """

    def __init__(self, table_name: str = "grace") -> None:
        """
        Initialize the CalcParams object with the given table name.

        Parameters
        ----------
        table_name : str, optional
            The name of the table where the parameters are stored. By default,
            it is set to "grace".
        """
        super(CalcParams, self).__init__(
            table_name=table_name,
        )
