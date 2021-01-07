import logging
import sys
from wmaee.core.common import Shell, LoggerMixin, OutputEmitter, EventHandler
from wmaee.utils import ase_to_pymatgen, to_pymatgen
from wmaee.potentials import construct_potcar
from wmaee.core.runner import vasp, get_vasp_configuration
from wmaee.core.types import *
from pymatgen.io.vasp import Poscar, Potcar, Kpoints, Incar, Vasprun
from pymatgen import Structure
from ase import Atoms
from uuid import uuid4
from os import getcwd
from io import TextIOWrapper
from glob import glob
from fnmatch import fnmatch
from os.path import exists, isfile, isdir, join, dirname
from tarfile import open as tar_open


def _write_input(structure: Union[Atoms, Collection[Atoms], Callable, Poscar], incar: Union[str, Incar],
                 kpoints: Optional[Union[None, str, Kpoints]] = None, potcar: Optional[Union[None, str, Potcar]] = None,
                 directory: Optional[Union[None, Directory]] = None, xc_func: Optional[str] = 'gga',
                 delete: Optional[bool] = False):
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
        logger.warning(
            'You did not specify a directory where you want to run VASP. I use "{}" for you.'.format(directory.name))

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
        elif isinstance(structure, Collection):
            logging.getLogger('wmaee.vasp._write_input').warning(
                'A list of structures was passed. I\'ll take the first one!')
            structure = next(iter(structure))
        elif isinstance(structure, Callable):
            logging.getLogger('wmaee.vasp._write_input').warning(
                'A list of structures was passed. I\'ll take the first one!')
            structure = structure(0)
        else:
            raise TypeError(type(structure))

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
            logger.warning(
                'No KPOINTS file was provided. That\'s not recommended. I\'ll proceed with a Gamma-Point only mesh')
            kpoints = Kpoints.gamma_automatic((1, 1, 1))

        if isinstance(kpoints, str):
            kpoints = Kpoints.from_file(kpoints)

        incar.write_file('INCAR')
        poscar.write_file('POSCAR')
        kpoints.write_file('KPOINTS')
        potcar.write_file('POTCAR')


class VASPInput(LoggerMixin):
    """
    A class to represent and write the input data needed by VASP.
    """

    def __init__(self, structure: Union[Atoms, Collection[Atoms]], incar: Union[str, Incar],
                 kpoints: Optional[Union[None, str, Kpoints]] = None, potcar: Optional[Union[None, str, Potcar]] = None,
                 xc_func: Optional[str] = 'gga'):
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

    def write(self, directory: Optional[Union[Directory, None]] = None):
        """
        Writes the input data to a certain directory
        :param directory: (str or working_directory) the directory where the input files are written to (default: None)
        """
        directory = directory if directory is not None else working_directory(getcwd())
        _write_input(self._structure, self._incar, kpoints=self._kpoints, potcar=self._potcar, xc_func=self._xc_func,
                     directory=directory)

    @staticmethod
    def from_directory(directory: Optional[Union[Directory, None]] = None):
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
                incar=Incar.from_file('INCAR'),
                structure=Poscar.from_file('POSCAR').structure,
                kpoints=Kpoints.from_file('KPOINTS'),
                potcar=Potcar.from_file('POTCAR')
            )


def write_input(inp, directory: Optional[Union[Directory, None]] = None):
    """
    Write a VASP input to a directory. Is equivalent to `inp.write(directory)`
    :param inp: (VASPInput) the VASPInput object
    :param directory: (str or working_directory) the directory where the input files are written to (default: None)
    """
    _write_input(inp.structure, inp.incar, kpoints=inp.kpoints, potcar=inp.potcar, xc_func=inp.xc_func,
                 directory=directory)


class VASPOutput(LoggerMixin):
    __create_key = str(uuid4())

    @classmethod
    def create(cls, vasprun: Vasprun):
        """
        Creates a instance of an VASPOutput object
        :param vasprun: (pymatgen.io.vasp.Vasprun) the parsed XML data of the run
        :return:
        """
        if not isinstance(vasprun, Vasprun):
            raise TypeError(type(vasprun))
        return VASPOutput(cls.__create_key, vasprun)

    def __init__(self, create_key: str, vasprun: Vasprun):
        """
        Internal method to create VASPOutput instances

        :param create_key: (str) the magic key of the class
        :param vasprun: (pymatgen.io.vasp.Vasprun) the parsed XML data
        """
        super(VASPOutput, self).__init__()
        assert (create_key == VASPOutput.__create_key), \
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


def _parse_archive(a : str, raise_exc: Optional[bool] = True) -> Union[VASPOutput, List[VASPOutput]]:
    """
    TODO: Add some documentations here
    :param a:
    :param raise_exc:
    :return:
    """
    def _zopen_fake(filename, mode='r'):
        return filename
    # we have to inject some custom function, to make pymatgen.io.vasp.Vasprun work also with file-like objects
    path_pattern = '*vasprun.xml'
    module_vasprun = sys.modules[Vasprun.__module__]
    module_potcar = sys.modules[Potcar.__module__]
    old_zopen = getattr(module_vasprun, 'zopen')
    old_get_potcars = getattr(Vasprun, 'get_potcars')

    with tar_open(a) as archive:
        vasp_run_candidates = [m for m in archive.getmembers() if fnmatch(m.name, path_pattern)]
        if len(vasp_run_candidates) < 1:
            raise FileNotFoundError('Could not find VASP output file "vasprun.xml"')

        try:
            out_handles = [archive.extractfile(cand) for cand in vasp_run_candidates]
            setattr(module_vasprun, 'zopen', _zopen_fake)
            setattr(module_potcar, 'zopen', _zopen_fake)

            out_data = []
            for out_h, cnd in zip(out_handles, vasp_run_candidates):
                potcar_member = archive.getmember(join(dirname(cnd.name), 'POTCAR'))
                potcar_handle = TextIOWrapper(archive.extractfile(potcar_member))
                potcar = Potcar.from_file(potcar_handle)

                def _get_potcars_fake(self, path):
                    return potcar

                setattr(Vasprun, 'get_potcars', _get_potcars_fake)
                out_data.append(VASPOutput.create(Vasprun(out_h)))
        except:
            if raise_exc:
                raise
        finally:
            setattr(module_vasprun, 'zopen', old_zopen)
            setattr(module_potcar, 'zopen', old_zopen)
            setattr(Vasprun, 'get_potcars', old_get_potcars)
    return out_data[0] if len(out_data) == 1 else out_data


def _parse_directory(d : working_directory, raise_exc = True) -> VASPOutput:
    """
    Parses a VASP calculation directory
    :param d: (working_directory)
    :param raise_exc: (bool)
    :return: (VASPOutput)
    """
    with d:
        if not exists('vasprun.xml'):
            raise FileNotFoundError('Could not find VASP output file "vasprun.xml"')
        try:
            vasprun = Vasprun('vasprun.xml')
        except Exception:
            if raise_exc:
                raise
        return VASPOutput.create(vasprun)


def parse_output(directory=None, raise_exc=True):
    """
    Parses the output from VASP. Searches for 'vasprun.xml' file
    :param directory: (str or working_directory) the directory where the output files are located (default: None)
    :return: (VASPOutput) the VASPOutput object
    """
    if directory is None:
        directory = working_directory(getcwd())

    if isinstance(directory, working_directory):
        return _parse_directory(directory, raise_exc=raise_exc)
    elif isinstance(directory, str):
        dirs = glob(directory)
        directories = [d for d in dirs if isdir(d)]
        files = [d for d in dirs if isfile(d)]

        data = []
        for d in directories:
            data.append((d, _parse_directory(working_directory(d), raise_exc=raise_exc)))
        for f in files:
            data.append((f, _parse_archive(f, raise_exc=raise_exc)))
        return data


def full_run(inp: VASPInput, directory: Optional[Union[Directory, None]] = None, cpus: Optional[int] = 2,
             show_output: Optional[bool] = True, application: Optional[Union[str, None]] = None,
             hostname: Optional[Union[str, None]] = None, partition: Optional[Union[str, None]] = None):
    """
    Run VASP in a given directory. The maximum number of cores is limited to six. If no directory is specified VASP
    will be executed on os.getcwd() directory
    :param inp: (VASPInput) the VASP input
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :param application: (str) the string of the application name (default: None)
    :param hostname: (str) the string of the current hostname (default: None)
    :param partition: (str) the partition name as a string (default: None)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    if isinstance(inp.structure, (Callable, list, tuple)):
        result = _vasp_interactive_internal(inp, directory=directory, cpus=cpus, show_output=show_output,
                                            application=application,
                                            hostname=hostname, partition=partition)

    else:
        write_input(inp, directory=directory)
        result = vasp(directory=directory, cpus=cpus, show_output=show_output, application=application,
                      hostname=hostname, partition=partition, return_stdout=False)

    result = result == 0
    if result:
        return parse_output(directory=directory)
    else:
        raise RuntimeError('VASP crashed')


def _vasp_interactive_internal(inp: VASPInput, directory: Optional[Union[Directory, None]] = None,
                               cpus: Optional[int] = 2, show_output: Optional[bool] = True,
                               return_stdout: Optional[bool] = False,
                               application: Optional[Union[str, None]] = None,
                               hostname: Optional[Union[str, None]] = None,
                               partition: Optional[Union[str, None]] = None,
                               callbacks: Optional[Collection[Callable]] = None,
                               output: Optional[Union[Collection[TextIO]]]=None):
    """
    Run VASP in a given directory. The maximum number of cores is limited to six. If no directory is specified VASP
    will be executed on os.getcwd() directory
    :param inp: (VASPInput) the VASP input
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :param application: (str) the string of the application name (default: None)
    :param hostname: (str) the string of the current hostname (default: None)
    :param partition: (str) the partition name as a string (default: None)
    :param callbacks: (list of callable) callback methods after each ionic step (default: None)
    :param output: (list of TextIO) list of streams where stdout will be forked to (default: None)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    logger = logging.getLogger('wmaee.core.vasp.VASPRunner')

    if directory is None:
        directory = working_directory(getcwd())

    if isinstance(directory, working_directory):
        pass
    elif isinstance(directory, str):
        directory = working_directory(directory)

    if isinstance(inp.structure, (list, tuple)):
        # it is not a collection thus therefore we do not have to do anything
        def mkgen():
            length = len(inp.structure)
            index = [0]
            iterator = iter(inp.structure)

            def gen_():
                is_last = index[0] + 1 == length
                index[0] += 1
                return is_last, next(iterator)

            return gen_

        generator = mkgen()
    elif isinstance(inp.structure, Callable):
        generator = inp.structure
    else:
        raise TypeError('For VASP Interactive "%s" is not allowed.' % type(inp.structure).__name__)
    stopcar = Incar(dict(LSTOP=True))
    # Check in INTERACTIVE = .True. in INCAR file
    inp.incar['INTERACTIVE'] = True
    inp.incar['ISYM'] = 0

    with directory:
        preamble, command, binary = get_vasp_configuration(application=application, hostname=hostname,
                                                           partition=partition)
        command = command.format(cores=cpus, binary=binary)
        change_command = 'cd {}'.format(getcwd())

        is_last, current_structure = generator()

        if is_last:
            inp.incar['INTERACTIVE'] = False
            inp.incar['ISYM'] = 0
            inp.incar['NSW'] = 0
            # Then we have to remove the interactive tag
        _write_input(structure=current_structure, incar=inp.incar, potcar=inp.potcar, kpoints=inp.kpoints,
                     xc_func=inp.xc_func)

        with Shell(stdout=sys.stdout if show_output else None, stderr=sys.stderr) as shell:
            shell.run(change_command)
            shell.run(preamble)
            position_callback = OutputEmitter('POSITIONS: reading from stdin')
            # create event handlers
            if callbacks is not None:
                handlers = [EventHandler('ionic_step_callback_handler_%i' % i, cb) for i, cb in enumerate(callbacks)]
                for handler in handlers:
                    position_callback.trigger.add_event_handler(handler)

            def set_new_positions():
                is_last, next_structure = generator()
                next_structure = to_pymatgen(next_structure)
                for atom in next_structure.frac_coords:
                    text = " ".join(map("{:19.16f}".format, atom))
                    shell._send_command(text+'\n', raw=True)
                if is_last:
                    stopcar.write_file('STOPCAR')
                    shell._send_command('exit\n', raw=True)
                    # Disable also this handle
                    position_callback.trigger -= set_new_positions

            # only register the handle if more than one structure is in our list
            position_callback.trigger += set_new_positions
            output_streams = [position_callback] if output is None else [position_callback] + list(output)
            output = shell.run(command, return_out=return_stdout, return_exit=True, output=output_streams)
            if callbacks is not None:
                for handler in handlers:
                    position_callback.trigger.remove_event_handler(handler)
        if return_stdout:
            _, exitcode, output = output
        else:
            _, exitcode = output

        if exitcode != 0:
            logger.warning('VASP did not execute successfully! Exit-Code: "{}"'.format(exitcode))

        return exitcode if not return_stdout else (exitcode, output)


def vasp_interactive(inp: VASPInput, directory: Optional[Union[Directory, None]] = None,
                     cpus: Optional[int] = 2, show_output: Optional[bool] = True,
                     return_stdout: Optional[bool] = False,
                     application: Optional[Union[str, None]] = None,
                     hostname: Optional[Union[str, None]] = None,
                     partition: Optional[Union[str, None]] = None,
                     callbacks: Optional[Collection[Callable]] = None,
                     output: Optional[Union[Collection[TextIO]]] = None):
    """
    Run VASP in a given directory. The maximum number of cores is limited to six. If no directory is specified VASP
    will be executed on os.getcwd() directory
    :param inp: (VASPInput) the VASP input
    :param directory: (str or working_directory) the directory where to execute vasp
    :param cpus: (int) the number of core used to run the VASP subprocess (default: 2)
    :param show_output: (bool) wether to print the output on sys.stdout, passed on to _read_output (default: True)
    :param return_stdout: (bool) wether to return the output as a list of strings. Passed on to _send_command (default: False)
    :param application: (str) the string of the application name (default: None)
    :param hostname: (str) the string of the current hostname (default: None)
    :param partition: (str) the partition name as a string (default: None)
    :param callbacks: (list of callable) callback methods after each ionic step (default: None)
    :param output: (list of TextIO) list of streams where stdout will be forked to (default: None)
    :return: (bool) or (bool, list of str) exitcode== 0 and output depending on the setting of propagate_stdout
    """
    return _vasp_interactive_internal(inp, directory=directory, cpus=cpus, show_output=show_output,
                                      return_stdout=return_stdout, application=application, hostname=hostname,
                                      partition=partition, callbacks=callbacks, output=output)
