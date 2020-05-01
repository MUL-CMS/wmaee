from wmaee.core.common import LoggerMixin
from wmaee.utils import to_pyiron
from typing import Union, Optional, Dict, Collection, NoReturn, List, Generator, Tuple
from ase import Atoms as AseAtoms
from pyiron.atomistics.structure.atoms import Atoms as IronAtoms
from pyiron.base.settings.generic import Settings
from pyiron.project import Project
from pymatgen import Structure, Lattice
from pandas import DataFrame
from uuid import uuid4
from os import getcwd
from numpy import ndarray
import logging


class PotentialException(Exception):
    pass


class LAMMPSCalculation(LoggerMixin):
    """
    A class used for creating LAMMPS jobs, and is built around pyiron
    """

    def __init__(self, structure: Union[Structure, AseAtoms, IronAtoms], potential: Union[str, None] = None,
                 name: Optional[Union[str, None]] = None, directory: Optional[Union[str, None]] = None):
        super(LAMMPSCalculation, self).__init__()
        self._structure = structure

        if directory is None:
            directory = 'lammps_%s' % uuid4().hex
            self.logger.info('No working directory was specified. I\'ve created "%s" for you' % directory)
        self._directory = directory

        if potential is not None:
            if not self.is_valid_potential(structure, potential):
                raise PotentialException
        self._potential = potential

        if name is None:
            name = 'lammps_calc_%s' % uuid4().hex
            self.logger.info('The calculation was no given a name hence I\'ve chosen "%s" for you' % name)

        self._name = name
        self._pyiron_settings_reference = Settings()
        self._modify_pyiron_config()
        self._project_handle = Project(self._directory)
        self._pyiron_job = self._project_handle.create_job(self._project_handle.job_type.Lammps, self._name)
        self._pyiron_job.structure = to_pyiron(self._structure)
        self._pyiron_job.potential = self._potential

    def _modify_pyiron_config(self):
        """
        Inserts the current directory into the pyiron configuration, so that pyiron can run jobs there
        """
        pwd = getcwd()
        if pwd not in self._pyiron_settings_reference._configuration['project_paths']:
            self._pyiron_settings_reference._configuration['project_paths'].append(pwd)

    @property
    def structure(self) -> Union[Structure, AseAtoms, IronAtoms]:
        return self._structure

    @structure.setter
    def structure(self, other: Union[Structure, AseAtoms, IronAtoms]) -> NoReturn:
        if not self.is_valid_potential(other, self._potential):
            self.logger.warning(
                'The potential "%s" is not valid for this structure. Consider to set a new potential' % self._potential)
        self._structure = other
        self._pyiron_job.structure = to_pyiron(self._structure)

    @property
    def potential(self) -> str:
        return self._potential

    @potential.setter
    def potential(self, other: str) -> NoReturn:
        if not self.is_valid_potential(self._structure, other):
            self.logger.warning(
                'The potential "%s" is not valid for this structure. Consider to set a new structure' % other)
        self._potential = other
        self._pyiron_job.potential = self._potential

    @classmethod
    def is_valid_potential(cls, structure: Union[Structure, AseAtoms, IronAtoms, None], potential: str) -> bool:
        available_potentials = list(cls.list_potentials(structure)['Name'])
        if potential not in available_potentials:
            logging.getLogger(cls.fullname()).warning('Potential "%s" is not a valid potential' % potential)
            return False
        else:
            return True

    @classmethod
    def list_potentials(cls, structure: Optional[Union[Structure, AseAtoms, IronAtoms, None]] = None,
                        fields: Optional[Collection[str]] = ('Name',)) -> DataFrame:
        """
        List all interatomic potentials for the current atomistic sturcture including all potential parameters.
        To quickly get only the names of the potentials youcan
        :param structure: (ase.Atoms, Structure or pyiron.Atoms) structure which should be calculated
        :param fields: (list of str) columns to select from the data
        :return: (pandas.DataFrame) DataFrame including all potential parameters.
        """
        from pyiron.lammps.potential import LammpsPotentialFile
        if structure is None:

            list_of_elements = []
        else:
            structure = to_pyiron(structure)
            list_of_elements = set(structure.get_chemical_symbols())

        list_of_potentials = LammpsPotentialFile().find(list_of_elements)
        if fields is not None:
            list_of_potentials = list_of_potentials[list(fields)]

        if list_of_potentials is not None:
            list_of_potentials = list_of_potentials.reset_index(drop=True)
            return list_of_potentials
        else:
            raise TypeError(
                "No potentials found for this kind of structure: ",
                str(list_of_elements),
            )

    def md(self, temperature: Optional[Union[None, float]] = None, pressure: Optional[Union[None, float]] = None,
           n_ionic_steps: Optional[int] = 1000, time_step: Optional[float] = 1.0, n_print: Optional[int] = 100,
           temperature_damping_timescale: Optional[int] = 100.0, pressure_damping_timescale: Optional[float] = 1000.0,
           seed: Optional[Union[int, None]] = None, tloop=None,
           initial_temperature: Optional[Union[float, None]] = None,
           langevin: Optional[bool] = False, delta_temp: Optional[Union[float, None]] = None,
           delta_press: Optional[Union[float, None]] = None) -> NoReturn:
        """
        Set an MD calculation within LAMMPS. Nosé Hoover is used by default.

        :param temperature: (None or float) Target temperature. If set to None, an NVE calculation is performed. It is required when the pressure is set or langevin is set
        :param pressure: (None or float) Target pressure. If set to None, an NVE or an NVT calculation is performed. (This tag will allow for a list in the future as it is done for calc_minimize())
        :param n_ionic_steps: (int) Number of ionic steps
        :param time_step: (float) Step size between two steps. In fs if units==metal
        :param n_print: (int)  Print frequency
        :param temperature_damping_timescale: (float) The time associated with the thermostat adjusting the temperature. (In fs. After rescaling to appropriate time units, is equivalent to Lammps' `Tdamp`.)
        :param pressure_damping_timescale: (float) The time associated with the barostat adjusting the temperature. (In fs. After rescaling to appropriate time units, is equivalent to Lammps' `Pdamp`.)
        :param seed: (int) Seed for the random number generation (required for the velocity creation)
        :param tloop: TODO: Find out what this parameter is good for
        :param initial_temperature: (None or float)  Initial temperature according to which the initial velocity field is created. If None, the initial temperature will be twice the target temperature (which would go immediately down to the target temperature as described in equipartition theorem). If 0, the velocity field is not initialized (in which case  the initial velocity given in structure will be used). If any other number is given, this value is going to be used for the initial temperature.
        :param langevin: (bool) Activate Langevin dynamics
        :param delta_temp: (float) Thermostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
        :param delta_press: (float) Barostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
        """
        self._pyiron_job.calc_md(temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps,
                                 time_step=time_step, n_print=n_print,
                                 temperature_damping_timescale=temperature_damping_timescale,
                                 pressure_damping_timescale=pressure_damping_timescale, seed=seed, tloop=tloop,
                                 initial_temperature=initial_temperature, langevin=langevin, delta_temp=delta_temp,
                                 delta_press=delta_press)

    def minimize(self, e_tol: Optional[float] = 0.0, f_tol: Optional[float] = 1e-4, max_iter: Optional[int] = 100000,
                 pressure: Optional[Union[Collection, float]] = None, n_print: Optional[int] = 100,
                 style: Optional[str] = 'cg') -> NoReturn:
        """
        Sets parameters required for minimization.
        :param: e_tol: (float) If the magnitude of difference between energies of two consecutive steps is lower than or equal to `e_tol`, the minimisation terminates. (Default is 0.0 eV.)
        :param f_tol: (float) If the magnitude of the global force vector at a step is lower than or equal to `f_tol`, the minimisation terminates. (Default is 1e-4 eV/angstrom.)  max_iter (int): Maximum number of minimisation steps to carry out. If the minimisation converges before `max_iter` steps, terminate at the converged step. If the minimisation does not converge up to `max_iter` steps, terminate at the `max_iter` step. (Default is 100000.)
        :param pressure: (float, list, tuple, numpy.ndarray) Pressure in GPa at which minimisation is to be carried out. If None, isochoric (constant volume) condition will be used. If a float, cell shape changes are allowed. in each of the three primary axes, but no shearing occurs. If array-like, (up to) the first six elements will be interpreted as the x, y, z, xy, xz, and yz components of the pressure tensor, respectively. If this component-wise mode is used, cell changes corresponding to any one element of the stress tensor can be selectively disabled by setting this element to None. (Default is None, run isochorically.)
        :param n_print: (int) Write (dump or print) to the output file every n steps (Default: 100)
        :param style: ('cg'/'sd'/other values from Lammps docs) The style of the numeric minimization, either conjugate gradient, steepest descent, or other keys permissible from the Lammps docs on 'min_style'. (Default is 'cg' -- conjugate gradient.)
        """
        self._pyiron_job.calc_minimize(e_tol=e_tol, f_tol=f_tol, max_iter=max_iter, pressure=pressure, n_print=n_print,
                                       style=style)

    def vcsgc(self, mu: Dict, ordered_element_list: List, target_concentration: Optional[Union[Dict, None]] = None,
              kappa: Optional[int] = 1000., mc_step_interval: Optional[int] = 100, swap_fraction: Optional[float] = 0.1,
              temperature_mc: Optional[Union[float, None]] = None, window_size: Optional[Union[float, None]] = None,
              window_moves: Optional[Union[int, None]] = None, temperature: Optional[Union[None, float]] = None,
              pressure: Optional[Union[None, float]] = None, n_ionic_steps: Optional[int] = 1000,
              time_step: Optional[float] = 1.0, n_print: Optional[int] = 100,
              temperature_damping_timescale: Optional[float] = 100.0,
              pressure_damping_timescale: Optional[float] = 1000.0, seed: Optional[Union[int, None]] = None,
              initial_temperature: Optional[Union[float, None]] = None, langevin: Optional[bool] = False) -> NoReturn:
        """
        Run variance-constrained semi-grand-canonical MD/MC for a binary system. In addition to VC-SGC arguments, all
        arguments for a regular MD calculation are also accepted.
        https://vcsgc-lammps.materialsmodeling.org

        Note:
            For easy visualization later (with `get_structure`), it is highly recommended that the initial structure
            contain at least one atom of each species.

        Warning:
            Assumes the units are metal, otherwise units for the constraints may be off.

        :param mu: (dict) A dictionary of chemical potentials, one for each element the potential treats, where the dictionary keys are just the chemical symbol. Note that only the *relative* chemical potentials are used here, such that the swap acceptance probability is influenced by the chemical potential difference between the two species (a more negative value increases the odds of swapping *to* that element.)
        :param ordered_element_list: (list) A list of the chemical species symbols in the order they appear in the definition of the potential in the Lammps' input file.
        :param target_concentration: (dict) A dictionary of target simulation domain concentrations for each species *in the potential*. Dictionary keys should be the chemical symbol of the corresponding species, and the sum of all concentrations must be 1. (Default is None, which runs regular semi-grand-canonical MD/MC without any variance constraint.)
        :param kappa: (int) Variance constraint for the MC. Larger value means a tighter adherence to the target concentrations. (Default is 1000.)
        :param mc_step_interval: (int) How many steps of MD between each set of MC moves. (Default is 100.) Must divide the number of ionic steps evenly.
        :param swap_fraction: (float) The fraction of atoms whose species is swapped at each MC phase. (Default is 0.1.)
        :param temperature_mc: (float) The temperature for accepting MC steps. (Default is None, which uses the MD temperature.)
        :param window_size: (float) The size of the sampling window for parallel calculations as a fraction of something unspecified in the VC-SGC docs, but it must lie between 0.5 and 1. (Default is None, window is determined automatically.)
        :param window_moves: (int) The number of times the sampling window is moved during one MC cycle. (Default is None, number of moves is determined automatically.)
        :param temperature: (None or float) Target temperature. If set to None, an NVE calculation is performed. It is required when the pressure is set or langevin is set
        :param pressure: (None or float) Target pressure. If set to None, an NVE or an NVT calculation is performed. (This tag will allow for a list in the future as it is done for calc_minimize())
        :param n_ionic_steps: (int) Number of ionic steps
        :param time_step: (float) Step size between two steps. In fs if units==metal
        :param n_print: (int)  Print frequency
        :param temperature_damping_timescale: (float) The time associated with the thermostat adjusting the temperature. (In fs. After rescaling to appropriate time units, is equivalent to Lammps' `Tdamp`.)
        :param pressure_damping_timescale: (float) The time associated with the barostat adjusting the temperature. (In fs. After rescaling to appropriate time units, is equivalent to Lammps' `Pdamp`.)
        :param seed: (int)  Seed for the random number generation (required for the velocity creation)
        :param initial_temperature: (None or float)  Initial temperature according to which the initial velocity field is created. If None, the initial temperature will be twice the target temperature (which would go immediately down to the target temperature as described in equipartition theorem). If 0, the velocity field is not initialized (in which case  the initial velocity given in structure will be used). If any other number is given, this value is going to be used for the initial temperature.
        :param langevin: (bool) (True or False) Activate Langevin dynamics
        """
        self._pyiron_job.calc_vcsgc(mu=mu, ordered_element_list=ordered_element_list,
                                    target_concentration=target_concentration, kappa=kappa,
                                    mc_step_interval=mc_step_interval, swap_fraction=swap_fraction,
                                    temperature_mc=temperature_mc, window_size=window_size, window_moves=window_moves,
                                    temperature=temperature, pressure=pressure, n_ionic_steps=n_ionic_steps,
                                    time_step=n_ionic_steps, n_print=n_print,
                                    temperature_damping_timescale=temperature_damping_timescale,
                                    pressure_damping_timescale=pressure_damping_timescale, seed=seed,
                                    initial_temperature=initial_temperature, langevin=langevin, job_name="")

    def static(self):
        self._pyiron_job.calc_static()

    def run(self):
        self._pyiron_job.run()

    @property
    def steps(self) -> ndarray:
        """
        Get's the number of step ids
        :return: (numpy.ndarray) the step ids
        """
        return self._pyiron_job['output/generic/time']

    @property
    def forces(self) -> Generator[Tuple[int, ndarray]]:
        """
        Yields for each step the step id and the corresponding forces
        :return: (Generator(tuple(int, numpy.ndarray)) the generator object
        """
        for step, fcs in zip(self.steps, self._pyiron_job['output/generic/forces']):
            yield step, fcs

    @property
    def final_forces(self):
        return self._pyiron_job['output/generic/time'][-1], self._pyiron_job['output/generic/forces'][-1]

    @property
    def positions(self) -> Generator[Tuple[int, ndarray]]:
        """
        Yields for each step the step id and the corresponding array of positions
        :return: (Generator(tuple(int, numpy.ndarray)) the generator object
        """
        for step, pos in zip(self.steps, self._pyiron_job['output/generic/positions']):
            yield step, pos

    @property
    def final_positions(self):
        return self._pyiron_job['output/generic/time'][-1], self._pyiron_job['output/generic/positions'][-1]

    @property
    def structures(self) -> Generator[Tuple[int, Structure]]:
        """
        Yields for each step the step id and the corresponding structure
        :return: (Generator(tuple(int, Structure)) the generator object
        """
        species_list = self._pyiron_job.structure.get_chemical_symbols().tolist()
        for step, cell, positions in zip(self.steps, self._pyiron_job['output/generic/cells'],  self._pyiron_job['output/generic/positions']):
            yield step, Structure(Lattice(cell), species_list, positions, coords_are_cartesian=True)

    @property
    def final_structure(self):
        species_list = self._pyiron_job.structure.get_chemical_symbols().tolist()
        return self._pyiron_job['output/generic/time'][-1], Structure(Lattice(self._pyiron_job['output/generic/cells'][-1]), species_list, self._pyiron_job['output/generic/positions'][-1], coords_are_cartesian=True)

    @property
    def volumes(self) -> Generator[Tuple[int, float]]:
        """
        Yields for each step the step id and the corresponding volumes
        :return: (Generator(tuple(int, float)) the generator object
        """
        for step, vol in zip(self.steps, self._pyiron_job['output/generic/volume']):
            yield step, vol

    @property
    def final_volume(self):
        return self._pyiron_job['output/generic/time'][-1], self._pyiron_job['output/generic/volume'][-1]

    @property
    def pressures(self) -> Generator[Tuple[int, ndarray]]:
        """
        Yields for each step the step id and the corresponding pressures
        :return: (Generator(tuple(int, numpy.ndarray)) the generator object
        """
        for step, pressure in zip(self.steps, self._pyiron_job['output/generic/pressures']):
            yield step, pressure

    @property
    def final_pressure(self) -> Tuple[int, ndarray]:
        """
        Get's the final pressure of the calculation
        :return: (tuple(int, numpy.ndarry)) the final step id and the pressure
        """
        return self._pyiron_job['output/generic/time'][-1], self._pyiron_job['output/generic/pressures'][-1]

    @property
    def energies(self) -> Generator[Tuple[int, float]]:
        """
        Yields for each step the step id and the corresponding energy
        :return: (Generator(tuple(int, float)) the generator object
        """
        for step, energy in zip(self.steps, self._pyiron_job['output/generic/energy_pot']):
            yield step, energy

    @property
    def final_energies(self) -> Tuple[int, float]:
        """
        Get's the final energy of the calculation
        :return: (tuple(int, float)) the last step number and corresponding energy
        """
        return self._pyiron_job['output/generic/time'][-1], self._pyiron_job['output/generic/energy_pot'][-1]