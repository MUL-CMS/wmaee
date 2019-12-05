import abc
import re
import logging
import tarfile
import json
from wmaee.extensions.common import LoggerMixin, remove_white
from os.path import isfile, join, isdir, exists
from os import listdir
from io import StringIO
from shutil import copyfileobj
from tempfile import NamedTemporaryFile

from pymatgen.io.vasp import Potcar, PotcarSingle

DEFAULT_CONFIG = 'defaults.json'
POTENTIAL_ARCHIVES = {}
DEFAULT_POTENTIALS = {}
def _get_configuration_directory():
    """
    Build the path the to configuration directory for this module
    :return: (str) the absolute path of the configuration directory
    """
    import os
    if 'WMAEE_CONFIG_DIR' not in os.environ:
        if exists('.config'):
            return join(os.getcwd(), '.config')
        else:
            raise RuntimeError('No configuration directory found')
    else:
        return os.environ['WMAEE_CONFIG_DIR']

class PotentialException(Exception):
    """
    A class to indicate that there is a Exception with POTCAR files
    """

    def __init__(self, msg):
        super(PotentialException, self).__init__(msg)



class PotentialArchive(LoggerMixin):
    """
    A class to work with VASP POTCAR potential archives
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, path, xc_func='gga'):
        self.path = path
        self.xc_func = xc_func

    @staticmethod
    def copy_file(fsrc, fdst, length=16 * 1024):
        """copy data from file-like object fsrc to file-like object fdst"""
        while 1:
            buf = fsrc.read(length)
            if not buf:
                break
            fdst.write(buf)

    @abc.abstractmethod
    def has_potential(self, identifier):
        """
        Returns a boolean wether the potential is available
        """

    @abc.abstractmethod
    def is_valid_potential(self, identifier):
        """ Returns a boolean wether the potential has a POTCAR and a PSCTR file """

    @abc.abstractmethod
    def check_valid_archive(self):
        """Returns a boolean wether the archive is valid (all potentials are valid and all default potentials are available """

    @abc.abstractmethod
    def potentials(self):
        """ Returns a list of all potential names """

    @abc.abstractmethod
    def potcar(self, identifier):
        """ Returns a file stream for the POTCAR file"""

    @abc.abstractmethod
    def psctr(self, identifier):
        """ Returns a file stream for the POTCAR file"""

    @abc.abstractmethod
    def get_potentials_for_element(self, element):
        """docstring"""

    def default_potential(self, element):
        """
        Returns the default POTCAR identifier as specified in "defaults.json"
        :param element: (str) the element a potential is needed for
        :return: (str) the identifier of the potential
        """
        global DEFAULT_POTENTIALS
        if element not in DEFAULT_POTENTIALS[self.xc_func]:
            raise PotentialException('No default potential configured for element "{}" for xc_type="{}"'
                                     .format(element, self.xc_func))
        else:
            return DEFAULT_POTENTIALS[self.xc_func][element]


class TarPotentialArchive(PotentialArchive):
    """
    Representes a POTCAR database packed in a tar.gz archive as obtained from the VASP webpage
    """

    def __init__(self, path):
        super(TarPotentialArchive, self).__init__(path)
        if not isfile(path):
            raise PotentialException('{1}: {0} is not a file'.format(path, self.__class__.__name__))
        else:
            if not path.split('.')[-1] in ['bz2', 'tar', 'gz', 'xz']:
                raise PotentialException('{1}: {0} is not a tar archive'.format(path, self.__class__.__name__))

        self._tarfile = tarfile.open(path, 'r:*')
        self._tarinfo = self._tarfile.getmembers()
        self._names = list(map(lambda info: info.name, self._tarinfo))
        self.check_valid_archive()
        #    raise PotentialException('{} is not a valid VASP potential archive'.format(path))

    def check_valid_archive(self):
        for potential in self.potentials():
            if not self.is_valid_potential(potential):
                self.logger.warning(
                    '{} default potential is corrupted or does not exist. Please do not use it!'.format(potential))
                return False

        default_potentials = DEFAULT_POTENTIALS[self.xc_func]

        for element, potential in default_potentials.items():
            if not self.is_valid_potential(potential):
                self.logger.warning(
                    '{} potential is corrupted or does not exist. Please do not use it!'.format(potential))
                return False
        return True

    def is_valid_potential(self, identifier):
        if self.has_potential(identifier):
            potcar_path = join(identifier, 'POTCAR')
            psctr_path = join(identifier, 'PSCTR')
            return potcar_path in self._names and psctr_path in self._names
        else:
            return False

    def has_potential(self, identifier):
        return identifier in self._names

    def potentials(self):
        return list(
            map(lambda info: info.name,
                list(filter(
                    lambda info: info.isdir(), self._tarinfo)
                )
                )
        )

    def get_potentials_for_element(self, element):
        return list(filter(lambda pot: pot.startswith(element), self.potentials()))

    def potcar(self, identifier):
        if self.is_valid_potential(identifier):
            potcar_path = join(identifier, 'POTCAR')
            if potcar_path in self._names:
                try:
                    member = list(filter(lambda inf: inf.name == potcar_path, self._tarinfo))[0]
                    file_obj = self._tarfile.extractfile(member)
                except:
                    raise PotentialException('An error occured while extracting {0}.'.format(potcar_path))
                else:
                    return file_obj
            else:
                raise PotentialException('The POTCAR file for the potential {0} was not found.'.format(identifier))
        else:
            return None

    def psctr(self, identifier):
        if self.is_valid_potential(identifier):
            psctr_path = join(identifier, 'PSCTR')
            if psctr_path in self._names:
                try:
                    member = list(filter(lambda inf: inf.name == psctr_path, self._tarinfo))[0]
                    file_obj = self._tarfile.extractfile(member)
                except:
                    raise PotentialException('An error occured while extracting {0}.'.format(psctr_path))
                else:
                    return file_obj
            else:
                raise PotentialException('The PSCTR file for the potential {0} was not found.'.format(identifier))
        else:
            return None


def _make_potential_archives():
    """
    Loads the default POTCAR table for all XC functionals configured
    """
    global POTENTIAL_ARCHIVES, DEFAULT_POTENTIALS
    default_potential_config = join(_get_configuration_directory(), DEFAULT_CONFIG)
    with open(default_potential_config, 'rb') as default_potential_config_file:
        default_potentials = json.load(default_potential_config_file)
        DEFAULT_POTENTIALS = default_potentials
    functionals = list(default_potentials.keys())
    resources_directory = _get_configuration_directory()
    found_directories = [f for f in listdir(resources_directory)
                         if f in functionals and isdir(join(resources_directory, f))]

    # Search for potential archives
    for functional_potential_directory in found_directories:
        # Search at first for .tar.gz files
        archives = [f for f in listdir(join(resources_directory, functional_potential_directory))
                    if f.endswith('.tar.gz')]
        archive_found = False
        for archive in archives:
            # Try to find a right potential archive
            try:
                functional_archive = TarPotentialArchive(
                    join(resources_directory, functional_potential_directory, archive))
            except PotentialException:
                continue
            else:
                # We found a valid potential archive
                POTENTIAL_ARCHIVES[functional_potential_directory] = functional_archive
                logging.getLogger().info('Found valid potential archive "{}"'.format(join(resources_directory,
                                                                                          functional_potential_directory,
                                                                                          archive)))
                archive_found = True
                break
        if not archive_found:
            try:
                functional_archive = DirectoryPotentialArchive(
                    join(resources_directory, functional_potential_directory))
            except PotentialException:
                logging.getLogger().warning(
                    'Could not find a potential archive for functional "{}"'.format(functional_potential_directory))
            else:
                POTENTIAL_ARCHIVES[functional_potential_directory] = functional_archive
                logging.getLogger().info('Found valid potential archive "{}"'.format(join(resources_directory,
                                                                                          functional_potential_directory)))


def _extract_species(poscar):
    """
    Extracts the correct order of the elements as Specified in the POSCAR to construct a POTCAR file
    :param poscar: (pymatge.io.vasp.Poscar) Poscar object from which to extract the elements names
    :return: (list of str) a list of the species names in the structure
    """
    from pymatgen.core.periodic_table import Element
    with StringIO(poscar.get_string()) as poscar:
        # Skip the first 5 lines
        for _ in range(5):
            poscar.readline()
        element_line = [remove_white(crumb) for crumb in poscar.readline().split(' ') if
                        remove_white(crumb) != '']
        try:
            for element in element_line:
                Element(element)
        except ValueError:
            return []
        else:
            return element_line


def potcar_from_string(string):
    """
    Parses a string and creates a pymatgen.io.vasp.Potcar object from it
    :param string: (str) the path to the POTCAR file
    :return: (pymagen.io.vasp.Potcar) the Potcar object
    """
    fdata = string

    potcar = Potcar()
    potcar_strings = re.compile(r"\n?(\s*.*?End of Dataset)",
                                re.S).findall(fdata)
    functionals = []
    for p in potcar_strings:
        single = PotcarSingle(p)
        potcar.append(single)
        functionals.append(single.functional)
    if len(set(functionals)) != 1:
        raise ValueError("File contains incompatible functionals!")
    else:
        potcar.functional = functionals[0]
    return potcar


def construct_potcar(poscar, xc_func='gga'):
    """

    :param poscar: (pymatgen.io.vasp.Poscar) the Poscar object for which the Potcar should be created
    :param xc_func: (str) the name of the xc functional (as configure) (default: 'gga')
    :return: (pymatgen.io.vasp.Potcar) Potcar instance representing the corresponding POTCAR file
    """
    if xc_func == 'pbe':
        xc_func = 'gga'
    archive = POTENTIAL_ARCHIVES[xc_func]
    functional = xc_func
    poscar_species = _extract_species(poscar)
    with NamedTemporaryFile() as final:
        for element in poscar_species:
            try:
                default_potential = archive.default_potential(element)
            except PotentialException:
                logging.getLogger().warning(
                    'No default POTCAR found for "{}" for element "{}"'.format(functional, element))
                default_potential = element
            copyfileobj(archive.potcar(default_potential), final)
        final.seek(0)
        potcar = potcar_from_string(final.read().decode('utf-8'))
        if not poscar_species == [p.element for p in potcar]:
            raise PotentialException('Something went wrong while constructing the POTCAR file')
    return potcar

_make_potential_archives()