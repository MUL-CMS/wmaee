# Forward utility functions and classes
__author__ = 'Dominik Gehringer'
from os.path import exists
from os import getcwd
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar, Vasprun
from uuid import uuid4
from subprocess import Popen, PIPE
from io import TextIOWrapper
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms
from pymatgen import Structure
from wmaee.extensions.common import working_directory, LoggerMixin, StringStream
from wmaee.extensions.potentials import construct_potcar
import logging
import shlex


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




def greet():
    """
    Greets the students, participating in the exercises
    """
    print('Herzlich willkommen zu den Übungen zu Werkstoffmodellierung auf atomarer Ebene (DFT-Teil). Es freut uns, dass du da bist!')
    print('Welcome to the exercises of the lecture Materials Modelling on Atomistic Scals. We\'re happy to have you here!')


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
            logger.warning('Not POTCAR was provided. I\'ll construct one for you')
            potcar = construct_potcar(poscar, xc_func=xc_func)
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
        return self.ionic_step[-1]['e_fr_energy']
    
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

