from wmaee.core.common import LoggerMixin
from wmaee.utils import to_pyiron
from typing import Union, Optional, List, Collection, NoReturn
from ase import Atoms as AseAtoms
from pyiron.atomistics.structure.atoms import Atoms as IronAtoms
from pyiron.base.settings.generic import Settings
from pyiron.project import Project
from pymatgen import Structure
from pandas import DataFrame
from uuid import uuid4
import logging


class PotentialException(Exception):
    pass


class LAMMPSInput(LoggerMixin):
    """
    A class used for creating LAMMPS jobs, and is built around pyiron
    """

    def __init__(self, structure: Union[Structure, AseAtoms, IronAtoms], potential: Union[str, None] = None,
                 name: Optional[Union[str, None]] = None, directory: Optional[Union[str, None]] = None):
        super(LAMMPSInput, self).__init__()
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
        self._project_handle = Project(self._directory)
        self._pyiron_job = self._project_handle.create_job(self._project_handle.job_type.Lammps, self._name)
        self._pyiron_job.structure = to_pyiron(self._structure)
        self._pyiron_job.potential = self._potential

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

    def md(self, temperature: Optional[None, float] = None, pressure: Optional[None, float] = None,
           n_ionic_steps: Optional[int] = 1000, time_step: Optional[float] = 1.0, n_print: Optional[int] = 100,
           temperature_damping_timescale: Optional[int] = 100.0, pressure_damping_timescale: Optional[float] = 1000.0,
           seed: Optional[int, None] = None, tloop=None, initial_temperature: Optional[float, None] = None,
           langevin: Optional[bool] = False, delta_temp: Optional[float, None] = None,
           delta_press: Optional[float, None] = None):
        """
        Set an MD calculation within LAMMPS. Nosé Hoover is used by default.

        :param temperature (None or float): Target temperature. If set to None, an NVE calculation is performed. It is required when the pressure is set or langevin is set
        :param pressure (None or float): Target pressure. If set to None, an NVE or an NVT calculation is performed. (This tag will allow for a list in the future as it is done for calc_minimize())
        :param n_ionic_steps (int): Number of ionic steps
        :param time_step (float): Step size between two steps. In fs if units==metal
        :param n_print (int):  Print frequency
        :param temperature_damping_timescale (float): The time associated with the thermostat adjusting the temperature. (In fs. After rescaling to appropriate time units, is equivalent to Lammps' `Tdamp`.)
        :param pressure_damping_timescale (float): The time associated with the barostat adjusting the temperature. (In fs. After rescaling to appropriate time units, is equivalent to Lammps' `Pdamp`.)
        :param seed (int):  Seed for the random number generation (required for the velocity creation)
        :param tloop:
        :param initial_temperature (None or float):  Initial temperature according to which the initial velocity field is created. If None, the initial temperature will be twice the target temperature (which would go immediately down to the target temperature as described in equipartition theorem). If 0, the velocity field is not initialized (in which case  the initial velocity given in structure will be used). If any other number is given, this value is going to be used for the initial temperature.
        :param langevin (bool): (True or False) Activate Langevin dynamics
        :param delta_temp (float): Thermostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
        :param delta_press (float): Barostat timescale, but in your Lammps time units, whatever those are. (DEPRECATED.)
        """
        pass
