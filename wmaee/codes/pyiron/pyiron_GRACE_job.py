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

# from tensorpotential.calculator import grace_fm
# from tensorpotential.calculator import TPCalculator


from ase.io import write, read
import io
import ast
from os.path import join
from ase import Atoms
from ase.constraints import ExpCellFilter
from ase.io.trajectory import Trajectory
from typing import Tuple

import pandas as pd
import numpy as np
from typing import Literal, Optional

from wmaee.units import EV_PER_ANG3_TO_GPA


class Grace(AtomisticGenericJob):
    """
    Custom job class for using Grace within the pyiron framework.

    This class provides methods to perform static calculations, structure relaxation,
    and molecular dynamics simulations using Grace UMLIPs.

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

        This property controls the maximum force tolerance (in eV/Å) used during ionic relaxation.

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
        Get the maximum number of iterations for relaxation or MD run.

        Returns
        -------
        int
            Maximum number of iterations.
        """
        return self._generic_input['max_iter']

    @max_iter.setter
    def max_iter(self, max_iter: int):
        """
        Set the maximum number of iterations for relaxation or MD run.

        Parameters
        ----------
        max_iter : int
            Maximum number of iterations.
        """
        self._generic_input['max_iter'] = max_iter


    ####################################################
    ### --- IMPLEMENTATION OF STATIC CALCULATION --- ###
    ####################################################
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
        grace_model : {"GRACE-1L-MP-r6", "GRACE-2L-MP-r5", "GRACE-2L-MP-r6",
                       "GRACE-FS-OAM", "GRACE-1L-OAM", "GRACE-2L-OAM",
                       "GRACE-FS-OMAT", "GRACE-1L-OMAT", "GRACE-2L-OMAT"}, optional
            The Grace model to use for the calculation. Default is "GRACE-2L-OMAT".
    
        Raises
        ------
        ValueError
            If the structure is not set before calling this method.
    
        Notes
        -----
        This method sets up the Grace calculator with the specified model and
        prepares the job for a static calculation. The structure must be set
        before calling this method.
    
        Returns
        -------
        None
        """
        if self.structure is None:
            raise ValueError(
                'Structure must be set before calling calc_static().')
        self._generic_input['calc_mode'] = 'static'
        self.input['model'] = grace_model
        
    def _write_calc_static(self):
        """
        Write the calculation script for performing static calculations.

        This method generates a Python script that uses Grace for static calculations
        and writes it to a file in the working directory.
        """
        model = self.input["model"]

        script = [
            'from tensorpotential.calculator import grace_fm',
            'from ase.io import read',
            'from ase.io.cif import write_cif',
            '',
            'initial_structure = read("structure.cif")',
            f'initial_structure.calc = grace_fm("{model}")',
            '',
            '# perform the static calculation',
            'energy_pot = initial_structure.get_total_energy()',
            'forces = initial_structure.get_forces()',
            'stresses = initial_structure.get_stress()',
            '',
            '# write the results to an output file',
            'with open("log.out", "w") as f:',
            '    f.write("step\tenergy_pot\tforces\tstresses\\n")',
            '    f.write(f"1\\t{energy_pot}\\t{forces.tolist()}\\t{stresses.tolist()}\\n")',
            '',
            '# write the final structure to a CIF file',
            'write_cif("final_structure.cif", initial_structure)',
        ]
        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

    def _parse_calc_static(
        self, 
        output: str = 'log.out', 
        skip: int = 0
    ) -> Optional[pd.DataFrame]:
        """
        Parse the output of a static calculation and return optimization results.
    
        This method reads a static calculation output file and extracts the 
        optimization steps and energy information, returning them as a DataFrame.
    
        Parameters
        ----------
        output : str, optional
            Name of the output file to parse. Default is 'log.out'.
        skip : int, optional
            Number of rows to skip at the start of the file before parsing. Default is 0.
    
        Returns
        -------
        pd.DataFrame or None
            A DataFrame containing the optimization results if parsing succeeds,
            otherwise None.
        """

        try:
            df = pd.read_csv(
                join(self.working_directory, output),
                sep=r"\t+", 
                engine='python'
            )
            # convert strings to np.array
            df['forces'] = df['forces'].apply(lambda x: np.array(ast.literal_eval(x)))
            df['stresses'] = df['stresses'].apply(lambda x: np.array(ast.literal_eval(x)))
            return df
        except:
            self.logger.warning(f'Cannot parse output from {output} file.')
            return None        

    
    ##########################################################
    ### --- IMPLEMENTATION OF MINIMIZATION CALCULATION --- ###
    ##########################################################

    def calc_minimize(
        self,
        grace_model: Literal[
            "GRACE-1L-MP-r6", "GRACE-2L-MP-r5", "GRACE-2L-MP-r6",
            "GRACE-FS-OAM", "GRACE-1L-OAM", "GRACE-2L-OAM",
            "GRACE-FS-OMAT", "GRACE-1L-OMAT", "GRACE-2L-OMAT"
        ] = "GRACE-2L-OMAT",
        algorithm: Literal["FIRE", "BFGS", "LBFGS"] = "FIRE",
        algorithm_kwargs: dict | None = None,
        relax_cell: bool = False,
        relax_cell_kwargs: dict | None = None,
        ionic_force_tolerance: float = 0.001,
        max_iter: int = 500,
        n_print : int = 1,
        save_path: str | None = "relax.traj",
    ) -> None:
        """
        Perform a structure relaxation (geometry optimization).
    
        Parameters
        ----------
        grace_model : {"GRACE-1L-MP-r6", "GRACE-2L-MP-r5", "GRACE-2L-MP-r6",
                       "GRACE-FS-OAM", "GRACE-1L-OAM", "GRACE-2L-OAM",
                       "GRACE-FS-OMAT", "GRACE-1L-OMAT", "GRACE-2L-OMAT"}, optional
            The GRACE model to use for the calculation. Default is "GRACE-2L-OMAT".
        algorithm : {"FIRE", "BFGS", "LBFGS"}, optional
            Ionic relaxation algorithm from ASE. Default is "FIRE".
        algorithm_kwargs : dict, optional
            Additional keyword arguments for the relaxation algorithm. Default is None.
        relax_cell : bool, optional
            If True, relax the simulation cell. Default is False.
        relax_cell_kwargs : dict, optional
            Additional keyword arguments for cell relaxation. Default is None.
        ionic_force_tolerance : float, optional
            Convergence criterion for ionic forces in eV/Å. Default is 0.001.
        max_iter : int, optional
            Maximum number of optimization steps. Default is 500.
        n_print : int, optional
            Frequency of logging output (structures/stresses). Default is 1.
        save_path : str or None, optional
            File path to save the relaxation trajectory. Default is "relax.traj".
    
        Returns
        -------
        None
        """

        
        super().calc_minimize(
            ionic_energy_tolerance=0,
            ionic_force_tolerance=ionic_force_tolerance,
            max_iter=max_iter,
        )

        # Saving the parameters in dictionary
        self._generic_input['calc_mode'] = 'minimize'
        self.input['model'] = grace_model
        self.input['algo'] = algorithm
        self.input['algo_kwargs'] = algorithm_kwargs or {}
        self.input['relax_cell'] = relax_cell
        self.input['relax_cell_kwargs'] = relax_cell_kwargs or {}
        self._generic_input['ionic_force_tolerance'] = ionic_force_tolerance
        self._generic_input['n_print'] = n_print
        self._generic_input['max_iter'] = max_iter
        self.input['save_path'] = save_path
    
    def _write_calc_minimize(self):
        """
        Write the calculation script for performing structure relaxation.

        This method generates a Python script that uses Grace for structure
        relaxation and writes it to a file in the working directory.
        """

        model = self.input['model'] # Grace model
        algo = self.input.get('algo', 'FIRE') # ASE optimizer
        if algo not in ['FIRE', 'BFGS', 'LBFGS']:
            raise ValueError(
                f"Unsupported optimizer: {algo}. Supported optimizers are: FIRE, BFGS, LBFGS.")
        # ASE optimizer kwargs
        algo_kwargs = self.input.get('algo_kwargs', {}).copy()
        save_path = self.input.get('save_path', 'relax.traj') # Add trajectory parameter
        algo_kwargs['trajectory'] = f"{save_path}"
        algo_kwargs['loginterval'] = self._generic_input['n_print']
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
        relax_cell_kwargs_str = ", ".join(
            f'{k}="{v}"' if isinstance(v, str) else f"{k}={v}"
            for k, v in relax_cell_kwargs.items()
        )

        # Set parameters for minimization
        params_run = f"fmax={self._generic_input['ionic_force_tolerance']}, steps={self._generic_input['max_iter']}" 

        # Save parameters in        
        self._generic_input['relax_cell'] = relax_cell
        self._generic_input['save_path'] = save_path
        
        # Somewhat hidden feature of adding custom fixes
        # TODO: Document? Make it a more "official" paramerer?
        fixes = self.input.get('fixes', [])
        if isinstance(fixes, str):
            fixes = [fixes]
        elif isinstance(fixes, list):
            pass
        else:
            self.logger.warning('Fixes are neither string or array of strings. Ignoring')
            fixes = []
        fix_imports = self.input.get('fix_imports', [])
        if isinstance(fix_imports, str):                
            fix_imports = [fix_imports]
        elif isinstance(fix_imports, list):
            pass
        else:
            self.logger.warning('Fix imports are neither string or array of strings. Ignoring')
            fix_imports = []
            
        script = [
            'from tensorpotential.calculator import grace_fm',
            'from ase.io import read',
            'from ase.io.cif import write_cif',
            'from ase.io.trajectory import Trajectory',
            f'from ase.optimize import {algo}',
        ]
        if relax_cell:
            script += ['from ase.constraints import ExpCellFilter as ECF']
        for f in fix_imports:
            script += [f]
        script += ['']

        script += [
            f'struct = read("structure.cif")',
            f'struct.calc = grace_fm("{model}")'
        ]
        for f in fixes:
            script.append(f'struct.set_constraint({f})')
        script += ['']
        
        if not relax_cell:
            script += [
                f'relaxation = {algo.upper()}(atoms=struct, {algo_kwargs_str}).run({params_run})',]
        else:
            script += [
                f'relaxation = {algo.upper()}(ECF(atoms=struct, {relax_cell_kwargs_str}), {algo_kwargs_str}).run({params_run})']
        script += ['']

        script += ['write_cif("final_structure.cif", struct)']

        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

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


    ################################################
    ### --- IMPLEMENTATION OF MD CALCULATION --- ###
    ################################################
    def calc_md(
        self,
        grace_model: Literal[
            "GRACE-1L-MP-r6", "GRACE-2L-MP-r5", "GRACE-2L-MP-r6",
            "GRACE-FS-OAM", "GRACE-1L-OAM", "GRACE-2L-OAM",
            "GRACE-FS-OMAT", "GRACE-1L-OMAT", "GRACE-2L-OMAT"
        ] = "GRACE-2L-OMAT",
        ensemble: Literal["NVT", "nvt", "Langevin", "langevin", "NPT", "npt"] = "NVT", 
        temperature: float = 300.0,
        starting_temperature: int | None = None,
        pressure: float | None = None,
        time_step: float = 2.0,
        n_ionic_steps: int = 100,
        n_print: int = 10,
        temperature_damping_timescale: float = 100.0,
        pressure_damping_timescale: float | None = None,
        trajectory: str = "md_run.traj",
        logfile: str = "md_run.log",
        random_seed: int | None = None,
        mask: Tuple[int] | np.ndarray | None = None
    ):
        """
        Run a molecular dynamics simulation using a GRACE foundation model.
    
        Parameters
        ----------
        grace_model : {"GRACE-1L-MP-r6", "GRACE-2L-MP-r5", "GRACE-2L-MP-r6",
                       "GRACE-FS-OAM", "GRACE-1L-OAM", "GRACE-2L-OAM",
                       "GRACE-FS-OMAT", "GRACE-1L-OMAT", "GRACE-2L-OMAT"}, optional
            Choice of GrACE foundation model (default is "GRACE-2L-OMAT").
    
        ensemble : {"NVT", "nvt", "Langevin", "langevin", "NPT", "npt"}, optional
            Type of statistical ensemble:
            - "NVT" or "Langevin": constant volume and temperature.
            - "NPT": constant pressure and temperature.
            Case-insensitive (default is "NVT").
    
        temperature : float, optional
            Target simulation temperature in Kelvin (default is 300.0).
    
        starting_temperature : int or None, optional
            Initial temperature of the system in Kelvin. If None, uses
            the target `temperature` (default is None).
    
        pressure : float or None, optional
            Target external pressure in GPa for NPT runs. Ignored if ensemble
            is NVT or Langevin (default is None).
    
        time_step : float, optional
            Integration time step in femtoseconds (default is 2.0 fs).
    
        n_ionic_steps : int, optional
            Number of MD steps to perform (default is 100).
    
        n_print : int, optional
            Frequency of log printing and trajectory writing (default is 10).
    
        temperature_damping_timescale : float, optional
            Thermostat coupling time in femtoseconds (default is 100.0 fs).
    
        pressure_damping_timescale : float or None, optional
            Barostat coupling time in femtoseconds. Only relevant for NPT runs
            (default is None).
    
        trajectory : str, optional
            Output trajectory file name (default is "md_run.traj").
    
        logfile : str, optional
            Output log file name (default is "md_run.log").
    
        random_seed : int or None, optional
            Seed for random velocity initialization. If None, system clock
            is used (default is None).
            
        mask : None, sequence of 3 ints (0/1), or (3,3) array-like, optional
            Specify which components of the computational box (strain tensor)
            are allowed to change under the barostat (default is None equalling
            to *No mask*). Only relevant for NPT simulations; it has no effect 
            for NVT / Langevin ensembles.
        """
        self._generic_input['calc_mode'] = 'md'
        ensemble = ensemble.lower()
        if ensemble in ['langevin', 'nvt']: #requested NVT, ignore NPT settings
            pressure = None
            mask = None
        if pressure == None:
            pressure_damping_timescale = None
        super().calc_md(
            temperature=temperature,
            pressure=pressure,
            n_ionic_steps=n_ionic_steps,
            time_step=time_step,
            n_print=n_print,
            temperature_damping_timescale=temperature_damping_timescale,
            pressure_damping_timescale=pressure_damping_timescale,
        )
        self.input['model'] = grace_model
        self.input['ensemble'] = ensemble.lower()
        self.input['trajectory'] = trajectory
        self.input['logfile'] = logfile 
        if starting_temperature is None:
            starting_temperature = temperature
        self.input['starting_temperature'] = starting_temperature
        self.input['random_seed'] = random_seed
        self.input['mask'] = mask

    def _write_calc_md(self):
        """
        Write the calculation script for performing molecular dynamics (MD).

        This method generates a Python script that uses CHGNet for running MD
        simulations and writes it to a file in the working directory.
        """
        # Somewhat hidden feature of adding custom fixes
        # TODO: Document? Make it a more "official" paramerer?
        fixes = self.input.get('fixes', [])
        if isinstance(fixes, str):
            fixes = [fixes]
        elif isinstance(fixes, list):
            pass
        else:
            self.logger.warning(
                'Fixes are neither string or array of strings. Ignoring')
            fixes = []
        fix_imports = self.input.get('fix_imports', [])
        if isinstance(fix_imports, str):
            fix_imports = [fix_imports]
        elif isinstance(fix_imports, list):
            pass
        else:
            self.logger.warning(
                'fix_imports are neither string or array of strings. Ignoring')
            fix_imports = []
        
        script = [
            'from tensorpotential.calculator import grace_fm',
            'from ase.md.velocitydistribution import MaxwellBoltzmannDistribution',
            'from ase.io import read',
            'from ase.io.cif import write_cif',
            'from ase.io.trajectory import Trajectory',
            'from ase import units'
        ]
        for f in fix_imports:
            script += [f]
        script += [
            '',
            'struct = read("structure.cif")',]
        for f in fixes:
            script.append(f'struct.set_constraint({f})')
        script += [
            f'calc = grace_fm("{self.input["model"]}")',
            'struct.calc = calc',
            '',
        ]
        if not self.input['random_seed'] == None:
            script += [
                'import numpy as np',
                f'np.random.seed({self.input["random_seed"]})',
            ]
        script += [
            f'MaxwellBoltzmannDistribution(struct, temperature_K={self.input["starting_temperature"]})',
            '',
        ]
        if self.input['ensemble'] == 'langevin':
            script += [
                'from ase.md.langevin import Langevin',
                f'dyn = Langevin(struct, timestep={self._generic_input["time_step"]}, temperature_K={self._generic_input["temperature"]}, friction={self._generic_input["temperature_damping_timescale"]})',
                ]
        else:
            pressure = self._generic_input["pressure"]
            if not pressure == None:
                pressure /= EV_PER_ANG3_TO_GPA
            pfactor = self._generic_input["pressure_damping_timescale"]
            if not pfactor == None:
                pfactor = f'{pfactor}*units.fs'
            script += [
                'from ase.md.npt import NPT',                
                f'dyn = NPT(struct, timestep={self._generic_input["time_step"]}*units.fs, temperature_K={self._generic_input["temperature"]}, ttime={self._generic_input["temperature_damping_timescale"]}*units.fs, externalstress={pressure}, pfactor={pfactor} , mask={self._generic_input["mask"]})',
                ]
        script += [
            f'traj = Trajectory("{self.input["trajectory"]}", "w", struct)',
            f'dyn.attach(traj.write, interval={self._generic_input["n_print"]})',
            '',
            f'logfile = open("{self.input["logfile"]}", "w")',
            'logfile.write("step\\tTime\\tEtot\\tEpot\\tEkin\\tT\\tVol\\ta\\tb\\tc\\tsigma_xx\\tsigma_yy\\tsigma_zz\\tsigma_yz\\tsigma_xz\\tsigma_xy\\n")',
            'def print_status():',
            '    values = [dyn.nsteps, dyn.nsteps*dyn.dt]',
            '    values += [struct.get_potential_energy()+struct.get_kinetic_energy()]',
            '    values += [struct.get_potential_energy(), struct.get_kinetic_energy()]',
            '    values += [struct.get_temperature()]',
            '    values += [struct.get_volume()]+list(struct.cell.lengths())',
            '    values += list(struct.get_stress())',
            '    status = "\\t".join(str(v) for v in values)',
            '    print(status)',
            '    logfile.write(status+"\\n")',
            '    logfile.flush()',
            f'dyn.attach(print_status, interval={self._generic_input["n_print"]})',
            '',
            f'dyn.run({self._generic_input["n_ionic_steps"]})',
            '',
            'write_cif("final_structure.cif", struct)'
        ]

        with open(join(self.working_directory, 'calc_script.py'), 'w') as f:
            f.writelines("\n".join(script))

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
            # df = pd.read_csv(
            #     join(self.working_directory, output),
            #     skiprows=1,
            #     sep=r'\s+'
            # )
            df = pd.read_csv(
                join(self.working_directory, output),
                sep=r"\t+", 
                engine='python'
            )
            return df
        except:
            self.logger.warning(f'Cannot parse output from {output} file.')
            return None

    ###############################
    ### --- COMMON FUNCTION --- ###
    ###############################
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
        if self._generic_input['calc_mode'] == 'static':
            self._write_calc_static()
        elif self._generic_input['calc_mode'] == 'minimize':
            self._write_calc_minimize()
        elif self._generic_input['calc_mode'] == 'md':
            self._write_calc_md()
        else:
            raise ValueError(
                "No calculation type set. Use calc_static(), calc_minimize() or calc_md() first.")
            
    def collect_output(self) -> None:
        """
        Collect and store the results of a calculation.
    
        This method parses the final structure, trajectory, and other output properties
        from the calculation and stores them in HDF5 format. The internal structure
        is updated, and various properties such as energies, structure, forces, and 
        stresses are saved to the output group.
    
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
        if self._generic_input['calc_mode'] in ['minimize', 'md']:
            # For minimize and md, we expect a trajectory file with the relaxed structure
            if self._generic_input['calc_mode'] == 'minimize':
                traj_file = self.input['save_path']
                traj = Trajectory(join(self.working_directory, traj_file))
            else:
                traj_file = self.input['trajectory']
                traj = read(join(self.working_directory, traj_file), index=":")
            try:
                with self.project_hdf5.open('output/generic') as h5out:
                    h5out['cells'] = np.array([s.get_cell() for s in traj])
                    h5out['positions'] = np.array([s.positions for s in traj])
                    h5out['stresses'] = np.array(
                        [voigt_to_tensor(s.get_stress()) for s in traj])
                    h5out['forces'] = np.array([s.get_forces() for s in traj])
            except:
                self.logger.warning(
                    f'Cannot parse output from {traj_file} file.')
        
        # Parse other properties
        if self._generic_input['calc_mode'] == 'static':
            df = self._parse_calc_static()  
            df['max_force'] = df['forces'].apply(lambda x: np.max(x))
            with self.project_hdf5.open('output/generic') as h5out:
                h5out['energy_pot'] = df['energy_pot'].values
                h5out['energy_tot'] = df['energy_pot'].values
                h5out['forces'] = df['forces'].values
                h5out['max_force'] = df['max_force'].values
                h5out['stresses'] = np.array(
                    [voigt_to_tensor(s)*EV_PER_ANG3_TO_GPA for s in df['stresses'].values])
                h5out['steps'] = df['step'].values
        elif self._generic_input['calc_mode'] == 'minimize':
            df = self._parse_calc_minimize()
            with self.project_hdf5.open('output/generic') as h5out:
                h5out['energy_pot'] = df['Energy'].values
                h5out['energy_tot'] = df['Energy'].values
                h5out['max_force'] = df['fmax'].values
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
