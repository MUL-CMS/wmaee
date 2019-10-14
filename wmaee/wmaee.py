# Forward utility functions and classes
__author__ = 'Dominik Gehringer'
from os.path import exists, join, isdir, isfile
from os import mkdir, getcwd, chdir, listdir
from shutil import rmtree, copyfileobj
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar, PotcarSingle, Vasprun
from uuid import uuid4
from subprocess import Popen, PIPE
from io import TextIOWrapper, StringIO
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen import Structure
from ase import Atoms
from tempfile import NamedTemporaryFile
from plotly.graph_objects import Figure, Scatter, layout
from pymatgen import Spin, Orbital, Structure, Lattice
import numpy as np
import abc
import json
import tarfile
import logging
import shlex
import re
from threading import Thread

DEFAULT_CONFIG = 'defaults.json'
EXERCISE_DIRECTORY='cms-exercise'
POTENTIAL_ARCHIVES = {}
DEFAULT_POTENTIALS = {}
__VASP_PREAMBLE = [
    'source /export/opt/intel/bin/compilervars.sh intel64',
    'ulimit -s unlimited',
    'FFTW_ROOT=/calc/dnoeger/software/lib/fftw-3.3.6',
    'MPICH_ROOT=/calc/dholec/Software/centos72/mpich-3.2-intel/instal',
    'VASP_ROOT=/calc/dholec/Software/centos72/vasp.5.4.1-pavel',
    'export PATH=${VASP_ROOT}-mpich-ifort-dynamic/bin:${MPICH_ROOT}/bin:${FFTW_ROOT}/bin:/calc/dnoeger/software/scph/src:${PATH}',
    'export LD_LIBRARY_PATH=${VASP_ROOT}-dynamic:${MPICH_ROOT}/lib:${FFTW_ROOT}/lib:${LD_LIBRARY_PATH}'
]


__VASP_COMMAND = 'mpirun -np {tasks} /calc/dholec/Software/centos72/vasp.5.4.1-pavel-mpich-ifort-dynamic/bin/vasp_std'
# Uncomment this section in order to configure the module for the use on "mul-hpc"
#__VASP_PREAMBLE = [
#    'module purge',
#    'module load intel',
#    'module load mvapich2/2.2',
#    'module load mkl',
#    'module load scalapack/2.0.2',
#    'ulimit -s unlimited'
#]
#__VASP_COMMAND = 'mpirun -np {tasks} /calc/dnoeger/software/vasp-intel-mvapich2-mkl/5.4.1/bin/vasp_std'
__END_MARK__ = '__VASP_FINISHED__'


def pymatgen_to_ase(structure):
    """
    Converts a `pymatgen.Structure` object to a `ase.Atoms` object
    :param structure: (pymatgen.Structure) the structure to convert
    :return: (ase.Atom) the converted object
    """
    return AseAtomsAdaptor.get_atoms(structure)

def ase_to_pymatgen(atoms):
    """
     Converts  a `ase.Atoms` object to a `pymatgen.Structure` object
    :param atoms: (ase.Atom) the structure to convert
    :return: (pymatgen.Structure) the converted object
    """
    return AseAtomsAdaptor.get_structure(atoms)

def remove_white(string):
    """
    Removes all whitespaces in a given string
    :param string: (str) the string
    :return: (str) a copy without whitespaces
    """
    whitespace = [' ', '\t', '\n']
    mystr = str(string)
    for removal in whitespace:
        mystr = mystr.replace(removal, '')
    return mystr

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


class StringStream(StringIO):
    """
    A class representing a dummy stream, which can be used to write data to and read from it
    """

    def __init__(self, string=''):
        """
        Create a string stream with a initial value
        :param string: (str) the initial value (default: "")
        """
        super(StringStream, self).__init__(initial_value=string)
        self._pos = 0
        self._remaining = 0
        self._length = 0

    def read(self, size=-1):
        """
        Performs a read operation on the StringStream object. Blocks if not enough data is available
        :param size: (int) number of characters to be read from the stream (default: -1)
        :return: (str) data read from the StringStream object
        """
        while self._remaining < size:
            sleep(0.01)
        result = super(StringStream, self).read()
        # Increase position, from current position seek( ..., 1)
        result_length = len(result)
        # Increase position, and cosume
        self._pos += result_length
        self._remaining = self._length - self._pos
        self.seek(self._pos)
        return result

    def write(self, s):
        """
        Write data to the StringStream object
        :param s: (str) the data
        """
        write_length = len(s)
        super(StringStream, self).write(s)
        # After write file is at the end
        # Seek from back and make it available
        self._length += write_length
        self._remaining = self._length - self._pos

    def readline(self, size=-1, block=True):
        """
        Reads a line from the StringStream object. Block if not a full line is available
        :param size: (int) number of characters to read (default: -1)
        :param block: (bool) wether to block until  a line is available (default: True)
        :return: (str) the data read from the StringStream object
        """
        if self.tell() != self._pos:
            self.seek(self._pos)
        result = super(StringStream, self).readline()
        result_length = len(result)

        self._pos += result_length
        self._remaining = self._length - self._pos
        # Seek new position
        self.seek(self._pos)
        # if block:
        #    while not result:
        #        result = super(StringStream, self).readline()
        #        sleep(0.025)
        return result


def _get_configuration_directory():
    """
    Build the path the to configuration directory for this module
    :return: (str) the absolute path of the configuration directory
    """
    import getpass
    return '/calc/{USER}/{DATA}/.config'.format(USER=getpass.getuser(), DATA=EXERCISE_DIRECTORY)

def greet():
    """
    Greets the students, participating in the exercises
    """
    print('Herzlich willkommen zu den Übungen zu Werkstoffmodellierung auf atomarer Ebene (DFT-Teil). Es freut uns, dass du da bist!')
    print('Welcome to the exercises of the lecture Materials Modelling on Atomistic Scals. We\'re happy to have you here!')

class LoggerMixin(object):
    """
    A mixin for logger
    """

    def __init__(self, **kwargs):
        super(LoggerMixin, self).__init__(**kwargs)

    @property
    def logger(self):
        return logging.getLogger(self.fullname())

    @classmethod
    def fullname(cls):
        name = '.'.join([
            cls.__module__,
            cls.__name__
        ])
        return name

    
class working_directory(object):

    """
    A convenience class which syntactic sugar, allowing the user to change the directories.
    Can also be nested.
    """

    def __init__(self, name=None, prefix=None, delete=False):
        """
        Constructs a working_directory object
        :param name: (str) name of the directory if None is given os.getcwd() will be used (default: None)
        :param prefix: (str) a prefix where to locate the directory (default: None)
        :param delete: (bool) wether to delete the directory after a with clause (default: False)
        """
        self._name = str(uuid4()) if not name else name
        self._delete = delete
        self._curr_dir = getcwd()
        self._active = False
        if prefix is not None:
            self._name = join(prefix, self._name)

    def __enter__(self):
        if not exists(self._name):
            mkdir(self._name)
        chdir(self._name)
        self._active = True

    def __exit__(self, exc_type, exc_val, exc_tb):
        chdir(self._curr_dir)
        if self._delete:
            rmtree(self._name)
        self._active = False

    @property
    def name(self):
        return self._name
    
    @property
    def active(self):
        return self._active
    
def view(structure, spacefill=True, show_cell=True, camera='perspective', particle_size=0.5, background='white', color_scheme='element', show_axes=True):
    """
    Constructs a nglview view to display a structure
    :param structure: (pymatgen.Structure, ase.Atoms) the structure to display
    :param spacefill: (bool) to set the atoms size to spacefilling (default: True)
    :param show_cell:  (bool) wether to draw the unit cell or not (default: True)
    :param camera: (str) which camera projections to use 'perspective' or 'orthographic' (default: 'perspective')
    :param particle_size: (float) the size of the atoms (default: 0.5)
    :param background: (str) the name of the background color (default: 'white')
    :param color_scheme: (str) the name of the coloring scheme. Please refer to nglview documentation (default: 'element')
    :param show_axes: (bool) wether to draw the coordinate system of the unit cell (default: True)
    :return: (nglview.View) view wrapper
    """
    try:
        import nglview
    except ImportError:
        raise ImportError('nglview is needed')
    if isinstance(structure, Atoms):
        atoms = structure
    elif isinstance(structure, Structure):
        atoms = pymatgen_to_ase(structure)
    else:
        raise TypeError
    view_  = nglview.show_ase(atoms)
    if spacefill:
        view_.add_spacefill(radius_type='vdw', color_scheme=color_scheme, radius=particle_size)
        # view.add_spacefill(radius=1.0)
        view_.remove_ball_and_stick()
    else:
        view_.add_ball_and_stick()
    if show_cell:
        if atoms.cell is not None:
            view_.add_unitcell()
    if show_axes:
        view_.shape.add_arrow([-2, -2, -2], [2, -2, -2], [1, 0, 0], 0.5)
        view_.shape.add_arrow([-2, -2, -2], [-2, 2, -2], [0, 1, 0], 0.5)
        view_.shape.add_arrow([-2, -2, -2], [-2, -2, 2], [0, 0, 1], 0.5)
    if camera != 'perspective' and camera != 'orthographic':
        print('Only perspective or orthographic is permitted')
        return None
    view_.camera = camera
    view_.background = background
    return view_

class PotentialException(Exception):
    """
    A class to indicate that there is a Exception with POTCAR files
    """
    def __init__(self, msg):
        super(PotentialException, self).__init__(msg)

        
class ThreadWithReturnValue(Thread):
    """
    Small wraper around threading.Thread which stores the return value of the executed function
    """

    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args,
                                        **self._kwargs)

    def join(self, *args, **kwargs):
        Thread.join(self, *args, **kwargs)
        return self._return


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
                self.logger.warning('{} default potential is corrupted or does not exist. Please do not use it!'.format(potential))
                return False

        default_potentials = DEFAULT_POTENTIALS[self.xc_func]

        for element, potential in default_potentials.items():
            if not self.is_valid_potential(potential):
                self.logger.warning('{} potential is corrupted or does not exist. Please do not use it!'.format(potential))
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
                functional_archive = DirectoryPotentialArchive(join(resources_directory, functional_potential_directory))
            except PotentialException:
                logging.getLogger().warning('Could not find a potential archive for functional "{}"'.format(functional_potential_directory))
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


def _construct_potcar(poscar, xc_func='gga'):
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
                logging.getLogger().warning('No default POTCAR found for "{}" for element "{}"'.format(functional, element))
                default_potential = element
            copyfileobj(archive.potcar(default_potential), final)
        final.seek(0)
        potcar = potcar_from_string(final.read().decode('utf-8'))
        if not poscar_species == [p.element for p in potcar]:
            raise PotentialException('Something went wrong while constructing the POTCAR file')
    return potcar


def _shell_alive(shell_handle):
    """
    Thest wether a shell is still open
    :param shell_handle: (subprocess.Popen) the process handle
    :return: (bool) a boolean flag indicating if the shell is still opened
    """
    return shell_handle and not isinstance(shell_handle.poll(), int)


def _send_command(shell, cmd, return_stdout=False, propagate_stdout=True):
    """
    Execute a command on a system shell and return the output
    :param shell: (subprocess.Popen, subprocess.PIPE, subprocess.PIPE, StringStream) the process handles and output pipes
    :param cmd: (str) the command to execute
    :param return_stdout: (bool) wether to return the commands output (default: False)
    :param propagate_stdout: (bool) wether to send the commands output to the shell_output pipe (default: True)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    shell_handle, shell_stdin, shell_stdout, shell_output = shell
    shell_input = shell_stdin
    if not _shell_alive(shell_handle):
        raise RuntimeError('Shell is not open')

    finish = '__local_command_finish_mark__'
    # self._shell_stdin.write(cmd + '\n')
    echo_cmd = 'echo {} $?\n'.format(finish)
    command = cmd.strip('\n')
    shell_input.write('; '.join([cmd, echo_cmd]))
    # Very important - flush() otherwise nothing will happen and the reader thread will get stuck in an infinite
    shell_input.flush()

    shout = []
    exit_status = 0
    for line in shell_stdout:
        if finish in line:
            try:
                mark, code = line.lstrip().rstrip().split(' ')
                exit_status = int(code)
                assert mark == finish
            except:
                shout.append(line)
                if propagate_stdout:
                    shell_output.write(line)
                    shell_stdout.flush()
            else:
                break
        else:
            shout.append(line)
            if propagate_stdout:
                shell_output.write(line)
                shell_stdout.flush()
    return exit_status == 0 if not return_stdout else (exit_status == 0, shout)


    
def _read_output(shell, log_file, show_output):
    """
    Reads the output of until __END_MARK__ is reached, terminates by returning the exit code
    :param shell: (subprocess.PIPE or StringStream) the shell output handle
    :param log_file: (file) the filehandle to the log file
    :param show_output: (bool) wether to print the output to sys.stdout
    :return: (int) the exitcode
    """
    _, _, _, f = shell
    line = f.readline()
    log_file.write(line)
    if show_output:
        print(line, end='')
    while __END_MARK__ not in line.strip():
        line = f.readline()
        if line:
            log_file.write(line)
            if show_output:
                print(line, end='')
        # If everything worked old_line should contain the exit status
    line, exitcode = line.strip().split(' ')
    log_file.flush()
    try:
        exitcode = int(exitcode)
    except ValueError:
        log_file.close()
        return None
    else:
        log_file.close()
        return exitcode


def _write_input(structure, incar, kpoints=None, potcar=None, directory=None, xc_func='gga', delete=False):
    """
    Writes input files for a VASP calculations
    :param structure: (ase.Atoms or pymatgen.Structure or str or pymatgen.io.vasp.Poscar) the structure to use in the calculation
    :param incar: (pymatgen.io.vasp.Incar or dict) the incar file
    :param kpoints: (pymatgen.io.vasp.Kpoints) the Kpoints instance. If None a Gamma Point k-mesh is used (default: None)
    :param potcar: (pymatgen.io.vasp.Potcar) the Potcar instance. If None a instance will be created automatically (default: None)
    :param directory: (str of working_directory) the directory where to run the calculations. If None os.getcwd() will be used (default: None)
    :param xc_func: (str) the name of the XC functional
    :param delete: (bool) will be passed on to the working_directory instance
    """
    logger = logging.getLogger('VASPRunner')
    if directory is None:
        directory = working_directory(getcwd())
        logger.warning('You did not specify a directory where you want to run VASP. I use "{}" for you.'.format(directory.name))
    
    if isinstance(directory, working_directory):
        pass
    elif isinstance(directory, str):
        directory = working_directory(directory, delete=delete)
        
    with directory:
        if isinstance(structure, Poscar):
            pass
        elif isinstance(structure, Atoms):
            structure = Poscar(ase_to_pymatgen(structure).get_sorted_structure())
        elif isinstance(structure, Structure):  
            structure = Poscar(structure.get_sorted_structure())
        elif isinstance(structure, str):
            structure = Poscar.from_file(structure)
        else:
            raise ValueError()
        poscar = structure
            
        incar = Incar(incar) if isinstance(incar, dict) else incar
        
        if potcar is None:
            logger.warn('Not POTCAR was provided. I\'ll construct one for you')
            potcar = _construct_potcar(poscar, xc_func=xc_func)
        elif isinstance(potcar, str):
            potcar = Potcar.from_file(potcar)
            
        if isinstance(incar, str):
            incar = Incar.from_file(incar)
        
        if kpoints is None:
            logger.warning('No KPOINTS file was provided. That\'s not recommended. I\'ll proceed with a Gamma-Point only mesh')
            kpoints = Kpoints.gamma_automatic((1,1,1))
            
        if isinstance(kpoints, str):
            kpoints = Kpoints.from_file(kpoints)
            
        incar.write_file('INCAR')
        poscar.write_file('POSCAR')
        kpoints.write_file('KPOINTS')
        potcar.write_file('POTCAR')
    

def _run_vasp_internal(directory=None, cpus=2, show_output=True, return_stdout=False):
    """
    Call VASP programm and obtain the output, by executing the contents of __VASP_PREAMBLE and __VASP_COMMAND in
    bash subprocess and piping the output
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    logger = logging.getLogger('VASPRunner')
    shell_handle = Popen(shlex.split('/bin/bash'), stdin=PIPE, stdout=PIPE)
    shell_stdin = TextIOWrapper(shell_handle.stdin, encoding='utf-8')
    shell_stdout = TextIOWrapper(shell_handle.stdout, encoding='utf-8')
    shell_output = StringStream()
    shell = (shell_handle, shell_stdin, shell_stdout, shell_output)
    
    if directory is None:
        directory = working_directory(getcwd())
    
    if isinstance(directory, working_directory):
        pass
    elif isinstance(directory, str):
        directory = working_directory(directory)
    
    with directory:
        log_file = '{}.vasp.log'.format(directory.name)
        log_file_fd = open(log_file, 'w')
        
        thr_read_log = ThreadWithReturnValue(target=_read_output, args=(shell, log_file_fd, show_output))
        thr_read_log.start()
        command = __VASP_COMMAND.format(tasks=cpus)
        echo_cmd = 'echo "{} $?"'.format(__END_MARK__)
        change_command = 'cd {}'.format(getcwd())
        _send_command(shell, change_command, return_stdout=False, propagate_stdout=False)
        for preamble in __VASP_PREAMBLE:
            _send_command(shell, preamble, return_stdout=False, propagate_stdout=False)
        
        output = _send_command(shell, command, return_stdout=return_stdout, propagate_stdout=True)
        _send_command(shell, echo_cmd, return_stdout=False, propagate_stdout=True)
        thr_read_log.join()
        
        if return_stdout:
            exitcode, output = output
        else:
            exitcode = output
        
        if not exitcode:
            logger.warning('VASP did not execute successfully! Exit-Code: "{}"'.format(exitcode))
            return
        
        return exitcode if not return_stdout else (exitcode, output)
        # Write input files
        
_make_potential_archives()


class VASPInput(LoggerMixin):
    """
    A class to represent and write the input data needed by VASP.
    """
    
    def __init__(self, incar, structure, kpoints, potcar=None, xc_func='gga'):
        """
        :param incar: (pymatgen.io.vasp.Incar or dict or str) the incar file
        :param structure: (ase.Atoms or pymatgen.Structure or str or pymatgen.io.vasp.Poscar) the structure to use in the calculation
        :param kpoints: (pymatgen.io.vasp.Kpoints or str) the Kpoints instance. If None a Gamma Point k-mesh is used (default: None)
        :param potcar: (pymatgen.io.vasp.Potcar or str) the Potcar instance. If None a instance will be created automatically (default: None)
        :param xc_func: (str) the name of the XC functional
        """
        super(VASPInput, self).__init__()
        self._incar = incar
        self._structure = structure
        self._kpoints = kpoints
        self._potcar = potcar
        self._xc_func = xc_func
    
    @property
    def incar(self):
        return self._incar
    
    @property
    def structure(self):
        return self._structure
    
    @property
    def kpoints(self):
        return self._kpoints
    
    @property
    def potcar(self):
        return self._potcar
    
    @property
    def xc_func(self):
        return self._xc_func
    
    @xc_func.setter
    def xc_func(self, value):
        if value not in ['lda', 'gga']:
            raise ValueError('Currently only "lda" and "gga" is allowed')
            
    @potcar.setter
    def potcar(self, value):
        self._potcar = value
        
    @incar.setter
    def incar(self, value):
        self._incar = value
        
    @kpoints.setter
    def kpoints(self, value):
        self._kpoints = value
        
    @structure.setter
    def structure(self, value):
        self._structure = value
        
    def write(self, directory=None):
        """
        Writes the input data to a certain directory
        :param directory: (str or working_directory) the directory where the input files are written to (default: None)
        """
        directory = directory if directory is not None else working_directory(getcwd())
        _write_input(self._structure, self._incar, kpoints=self._kpoints, potcar=self._potcar, xc_func=self._xc_func, directory=directory)
        
    @staticmethod
    def from_directory(directory=None):
        """
        Parses the input data from a given directory. Useful to use a previous calculations
        :param directory: (str or working_directory) the directory where the input files are located (default: None)
        :return: (VASPInput) the object containing the parsed data
        """
        if directory is None:
            directory = getcwd()
        directory = directory if isinstance(directory, working_directory) else working_directory(directory)
        with directory:
            files = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR']
            for f in files:
                if not exists(f):
                    raise FileNotFoundError('Could not find VASP input file "{}"'.format(f))
            return VASPInput(
                Incar.from_file('INCAR'),
                Poscar.from_file('POSCAR').structure,
                Kpoints.from_file('KPOINTS'),
                Potcar.from_file('POTCAR')
            )
    
def write_input(inp, directory=None):
    """
    Write a VASP input to a directory. Is equivalent to `inp.write(directory)`
    :param inp: (VASPInput) the VASPInput object
    :param directory: (str or working_directory) the directory where the input files are written to (default: None)
    """
    _write_input(inp.structure, inp.incar, kpoints=inp.kpoints, potcar=inp.potcar, xc_func=inp.xc_func, directory=directory)

def vasp(directory=None, cpus=2, show_output=True, return_stdout=False):
    """
    Run VASP in a given directory. The maximum number of cores is limited to six. If no directory is specified VASP
    will be executed on os.getcwd() directory
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    if cpus < 0:
        cpus = 1
    elif cpus > 6:
        cpus = 6
        logging.getLogger('VASPRunner').warning('Since you have to share our nodew with you colleagues we limited the maximum number of core to six!\n'
                                                'Weil Ihr euch die Systemrssourcen teilen müsst, haben wir die Anzahl der maximalen Kerne pro Job auf sechs limitiert!')
    else:
        cpus = cpus
    return _run_vasp_internal(directory=directory, cpus=cpus, show_output=show_output, return_stdout=return_stdout)


class VASPOutput(LoggerMixin):

    __create_key = str(uuid4())

    @classmethod
    def create(cls, vasprun):
        """
        Creates a instance of an VASPOutput object
        :param vasprun: (pymatgen.io.vasp.Vasprun) the parsed XML data of the run
        :return:
        """
        if not isinstance(vasprun, Vasprun):
            raise TypeError(type(vasprun))
        return VASPOutput(cls.__create_key, vasprun)

    def __init__(self, create_key, vasprun):
        assert(create_key == VASPOutput.__create_key), \
            "VASPOutput objects must be created using VASPOutput.create"
        self._vasprun = vasprun
        
    @property
    def parameters(self):
        return self._vasprun.parameters
    
    @property
    def structures(self):
        return self._vasprun.structures
    
    @property
    def initial_structure(self):
        return self._vasprun.initial_structure
    
    @property
    def final_structure(self):
        return self._vasprun.final_structure
    
    @property
    def fermi_energy(self):
        return self._vasprun.efermi
    
    @property
    def ionic_steps(self):
        return self._vasprun.ionic_steps
    
    @property
    def final_stress(self):
        return self.ionic_steps[-1]['stress']
    
    @property
    def stresses(self):
        for step in self.ionic_steps:
            yield step['stress']
            
    @property
    def final_forces(self):
        return self.ionic_steps[-1]['forces']
    
    @property
    def forces(self):
        for step in self.ionic_steps:
            yield step['forces']
            
    @property
    def final_energy(self):
        return self.ionic_steps[-1]['e_wo_entrp']
    
    @property
    def energies(self):
        for step in self.ionic_steps:
            yield step['e_wo_entrp']
            
    @property
    def final_electronic_entropy(self):
        return self.ionic_steps[-1]['e_0_energy']
    
    @property
    def electronic_entropies(self):
        for step in self.ionic_steps:
            yield step['e_0_energy']
            
    @property
    def final_free_energy(self):
        return self.ionic_steps['e_fr_energy']
    
    @property
    def free_energies(self):
        for step in self.ionic_steps:
            yield step['e_fr_energy']
            
    @property
    def total_dos(self):
        return self._vasprun.tdos
    
    @property
    def partial_dos(self):
        return self._vasprun.pdos
    
    @property
    def eigenvalues(self):
        return self._vasprun.eigenvalues

    @staticmethod
    def from_directory(directory=None):
        return parse_output(directory=directory)
    

def parse_output(directory=None):
    """
    Parses the output from VASP. Searches for 'vasprun.xml' file
    :param directory: (str or working_directory) the directory where the output files are located (default: None)
    :return: (VASPOutput) the VASPOutput object
    """
    if directory is None:
        directory = working_directory(getcwd())
    
    if isinstance(directory, working_directory):
        pass
    elif isinstance(directory, str):
        directory = working_directory(directory)
        
    with directory:
        if not exists('vasprun.xml'):
            raise FileNotFoundError('Could not find VASP output file "vasprun.xml"')
        vasprun = Vasprun('vasprun.xml')
        return VASPOutput.create(vasprun)

def full_run(inp, directory=None, cpus=2, show_output=True):
    write_input(inp, directory=directory)
    result = vasp(directory=directory, cpus=cpus, show_output=show_output)
    if result:
        return parse_output(directory=directory)
    else:
        raise RuntimeError('VASP crashed')

def plot_total_dos(output, efermi=False, erange=None):
    """
    Plot the total density of states from a given output.
    :param output: (VASPOutput) the instance of VASPOutput
    :param efermi: (bool) weth to rescale the x-axis such that the Fermi level is zero (default: False)
    :param erange: (float, float) the energy range to clip the plot. If None the whole range will be used (default: None)
    :return: (plotly.graph_objects.Figure) the plot
    """
    # Non spin polarized case
    dos = output.total_dos
    fig = Figure()
    spins = (Spin.up, Spin.down)
    max_val = []
    min_val = []
    for i, spin in enumerate(spins):
        density = dos.densities[spin]*spin.value
        max_val.append(np.amax(density))
        min_val.append(np.amin(density))
        energies = dos.energies
        energies = energies if not efermi else energies - dos.efermi
        scatter = Scatter(x=energies, y=density, name=spin.name)
        fig.add_trace(scatter)
        if i >= len(dos.densities) - 1:
            break
            
    fermi_energy = dos.efermi if not efermi else 0.0
    min_val, max_val = min(min_val), max(max_val)
    fig.update_layout(shapes=[layout.Shape(type="line", x0=fermi_energy, y0=min_val, x1=fermi_energy,y1=max_val)])
    if erange is not None:
        fig.update_xaxes(range=erange)
    return fig

def plot_projected_dos(output, combine_spins=False, efermi=False, orbitals=None, sum_density=True, combine_orbitals=None, erange=None):
    """
    Plot the projected density of states
    :param output: (VASPOutput) the VASP output object
    :param combine_spins: (bool) wether to sum up the up and down spin density for spin polarized calculations (default: False)
    :param efermi: (bool) wether to rescale the x-axis such that the Fermi level is zero (default: False)
    :param orbitals: (list or tuple of str) which orbitals to display e.g ('s', 'p', 'd') or ('px', 'py'). None means all (default: None)
    :param sum_density: (bool) wether to sum up the densities for the individual orbital
    :param combine_orbitals:
    :param erange: (float, float) the energy range to clip the plot. If None the whole range will be used (default: None)
    :return: (plotly.graph_objects.Figure) the plot
    """
    pdos = output.partial_dos[-1]
    energies = output.total_dos.energies
    energies = energies if not efermi else energies - output.total_dos.efermi
    fig = Figure()
    max_val = []
    min_val = []
    sd = []
    
    figure_data = []
    flatten = lambda l: [item for sublist in l for item in sublist]
    for orbital, data in pdos.items():
        if orbitals is not None:
            orbital_str = [o.name if isinstance(o, Orbital) else o for o in orbitals]
            if not any([orbital.name.startswith(orb_str) for orb_str in orbital_str ]):
                continue
                
        density = [(Spin.up, np.sum(list(data.values()), axis=0))] if combine_spins else [(spin, spin.value*den) for spin, den in data.items()]
        max_val.append(np.amax([d for _, d in density]))
        min_val.append(np.amin([d for _, d in density]))
        
        if combine_spins:
            figure_data.append(dict(x=energies, y=density[0][1], name=orbital.name, orbital=orbital, spin=Spin.up))
            sd.append((Spin.up, density[0][1]))
        else:
            for s,d in density:
                figure_data.append(dict(x=energies, y=d, name='{}-{}'.format(orbital.name, s.name), orbital=orbital, spin=s))
                sd.append((s,d))
    if combine_orbitals is not None:
        combinations = []
        
        for combination in combine_orbitals:
            current_combination = []
            for dict_ in figure_data:
                orbital_str = [o.name if isinstance(o, Orbital) else o for o in combination]
                orbital = dict_['orbital']
                if any([orbital.name.startswith(orb_str) for orb_str in orbital_str ]):
                    current_combination.append(dict_)
            combinations.append(current_combination)
        # summ up all with the same sping
        for dat in flatten(combinations):
            for i, f in enumerate(figure_data):
                if dat['orbital'] == f['orbital'] and dat['spin'] == f['spin']:
                    figure_data.pop(i)
                    break
        combined_data = []
        for combination in combinations:
            sum_up, sum_down = [], []
            name_up, name_down = '', ''
            xaxis = None
            for orbital in [o for o in combination if o['spin']==Spin.down]:
                sum_down.append(orbital['y'])
                name_down += orbital['orbital'].name + '-'
                xaxis = orbital['x']
            for orbital in [o for o in combination if o['spin']==Spin.up]:
                sum_up.append(orbital['y'])
                name_up += orbital['orbital'].name + '-'
                xaxis = orbital['x']
            name_up += 'up'
            name_down += 'down'
            # If combine_spins only Spin.up is used
            if len(sum_up) > 0:
                sum_up = np.sum(sum_up, axis=0)
                combined_data.append(dict(
                    x=xaxis,
                    name=name_up,
                    y=sum_up
                ))
            
            if not combine_spins:
                # Check if it was spin polarized calulation at all
                if len(sum_down) > 0:
                    sum_down = np.sum(sum_down, axis=0)
                    combined_data.append(dict(
                        x=xaxis,
                        name=name_down,
                        y=sum_down
                    ))
    else:
        combined_data = []
        
    all_data = combined_data + figure_data
    for d in all_data:
        fig.add_trace(Scatter(x=d['x'], y=d['y'], name=d['name']))
        max_val.append(np.amax(d['y']))
        min_val.append(np.amin(d['y']))
    if sum_density:
        if combine_spins:
            sd_up = np.sum([d for s,d in sd if s == Spin.up], axis=0)
            max_val.append(np.amax(sd_up))
            min_val.append(np.amin(sd_up))
            fig.add_trace(Scatter(x=energies, y=sd_up, name='sum'))
        else:
            sd_up = np.sum([d for s,d in sd if s == Spin.up], axis=0)
            max_val.append(np.amax(sd_up))
            min_val.append(np.amin(sd_up))
            fig.add_trace(Scatter(x=energies, y=sd_up, name='sum-up'))
            if len([d for s,d in sd if s == Spin.down]) > 0:
                sd_down = np.sum([d for s,d in sd if s == Spin.down], axis=0)
                max_val.append(np.amax(sd_down))
                min_val.append(np.amin(sd_down))
                fig.add_trace(Scatter(x=energies, y=sd_down, name='sum-down'))
    fermi_energy = output.total_dos.efermi if not efermi else 0.0
    min_val, max_val = min(min_val), max(max_val)
    fig.update_layout(shapes=[layout.Shape(type="line", x0=fermi_energy, y0=min_val, x1=fermi_energy,y1=max_val)])
    if erange is not None:
        fig.update_xaxes(range=erange)
    return fig
