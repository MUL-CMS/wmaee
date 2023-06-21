
import io
import os
import re
import sys
import time
import subprocess
import contextlib
import numpy as np
import ase.io.abinit
from frozendict import frozendict
from wmaee.core.utils import merge
from wmaee.units import HARTREE_TO_EV
from ase.calculators.calculator import FileIOCalculator, CalculationFailed, all_changes, Calculator, CalculatorSetupError, kpts2mp
from typing import Dict, Tuple, Any, TextIO, Union

ABINIT_PSEUDO_REGEX = re.compile(r"\"?(.*)[\",]?")
FLOAT_REGEX = r"[+-]?\d*\.\d+[eE]?[+-]?\d*"
INTEGER_REGEX = r"[+-]?\d+"
ABINIT_ENERGY_UNITS = '|'.join(("eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs"))


def identity(x):
    return x


def times(val: float):
    return lambda v: val*v


def parse_pseudo(path: str) -> Dict[str, str]:
    path, *_ = path.split(',')
    *_, xc, pps = os.path.basename(path).split('.')
    return dict(xc=xc, pps=pps)


def parse_ngkpt(*nums) -> Dict[str, Tuple[int, ...]]:
    return dict(ngkpt=tuple(map(int, nums)))


def parse_shiftk(*nums) -> Dict[str, Tuple[float, ...]]:
    return dict(shiftk=tuple(map(float, nums)))


from_hartree = times(HARTREE_TO_EV)
from_rydberg = times(HARTREE_TO_EV/2)
ENERGY_CONVERSION = dict(ev=identity, hartrees=from_hartree, hartree=from_hartree, ha=from_hartree, rydbergs=from_rydberg, rydberg=from_rydberg, ry=from_rydberg)

ABINIT_ENERGY_SCALAR_REGEX = re.compile(rf"({FLOAT_REGEX}|{INTEGER_REGEX})\s*({ABINIT_ENERGY_UNITS})?")


def make_parser(key: str, conversion, default: str):
    def parse_value(value: str, unit: str = default) -> Dict[str, float]:
        unit = unit.lower() if unit is not None else None

        try:
            value = float(value)
        except:
            raise ValueError(value)
        return {key: conversion.get(unit)(value) if unit is not None else value}

    return parse_value


ABINIT_KEYS = frozendict(
    pseudos=(parse_pseudo, ABINIT_PSEUDO_REGEX),
    ecut=(make_parser("ecut", ENERGY_CONVERSION, "eV"), ABINIT_ENERGY_SCALAR_REGEX),
    toldfe=(make_parser("toldfe", ENERGY_CONVERSION, "eV"), ABINIT_ENERGY_SCALAR_REGEX),
    nshiftk=(lambda v: dict(nshiftk=int(v.strip())), re.compile(f"({INTEGER_REGEX})")),
    kptopt=(lambda v: dict(kptopt=int(v.strip())), re.compile(f"({INTEGER_REGEX})")),
    nsppol=(lambda v: dict(nsppol=int(v.strip())), re.compile(f"({INTEGER_REGEX})")),
    ixc=(lambda v: dict(ixc=int(v.strip())), re.compile(f"({INTEGER_REGEX})")),
    ngkpt=(parse_ngkpt, re.compile(rf"({INTEGER_REGEX})\s+({INTEGER_REGEX})\s+({INTEGER_REGEX})")),
    shiftk=(parse_shiftk, re.compile(rf"({FLOAT_REGEX})\s+({FLOAT_REGEX})\s+({FLOAT_REGEX})"))
)


def parse_key(key: str, value: str):
    if key in ABINIT_KEYS:
        parser, regex = ABINIT_KEYS.get(key)
        return parser(*regex.match(value).groups())


def parse_file(handle):
    buf = io.StringIO()
    curr_key = None
    is_multiline_regex = re.compile(r"^([a-zA-Z]+)\s*$")
    is_command_regex = re.compile(r"^([a-zA-Z]+)")
    for line in handle:
        if line.startswith('#'):
            continue
        is_command = is_command_regex.match(line)
        is_multiline_command = is_multiline_regex.match(line)

        if is_command is None:
            buf.write(line)
        else:
            # check if there is something in the buffer
            if buf.tell() and buf.getvalue().strip():
                result = parse_key(curr_key, buf.getvalue())
                buf = io.StringIO()
                if result is not None:
                    yield result
            curr_key = is_command.group(0)
            if not is_multiline_command:
                _, value = line.split(' ', maxsplit=1)
                if '#' in value:  # strip off comments
                    value, _ = value.split('#', maxsplit=1)
                result = parse_key(curr_key, value)
                if result is not None:
                    yield result


def parse_abinit_input(path: Union[str, TextIO]) -> Dict[str, Any]:
    handle = open(path) if isinstance(path, str) else contextlib.nullcontext(path)
    with handle:
        return merge(*parse_file(handle))


def get_abinit_version(command):
    txt = subprocess.check_output([command, '--version']).decode('ascii')
    # This allows trailing stuff like betas, rc and so
    m = re.match(r'\s*(\d\.\d\.\d)', txt)
    if m is None:
        raise RuntimeError('Cannot recognize abinit version. '
                           'Start of output: {}'
                           .format(txt[:40]))
    return m.group(1)


class Abinit(FileIOCalculator):
    """Class for doing ABINIT calculations.

    The default parameters are very close to those that the ABINIT
    Fortran code would use.  These are the exceptions::

      calc = Abinit(label='abinit', xc='LDA', ecut=400, toldfe=1e-5)
    """

    implemented_properties = ['energy', 'forces', 'stress', 'magmom']
    ignored_changes = {'pbc'}  # In abinit, pbc is always effectively True.
    command = 'abinit < PREFIX.files > PREFIX.log'
    discard_results_on_any_change = True

    default_parameters = dict(
        xc='LDA',
        smearing=None,
        kpts=None,
        raw=None,
        pps='fhi')

    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='abinit', atoms=None, pp_paths=None,
                 v8_legacy_format=None, output='-',
                 **kwargs):
        """Construct ABINIT-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'abinit'.

        Examples
        ========
        Use default values:

        >>> h = Atoms('H', calculator=Abinit(ecut=200, toldfe=0.001))
        >>> h.center(vacuum=3.0)
        >>> e = h.get_potential_energy()

        """

        self.v8_legacy_format = v8_legacy_format
        self.pp_paths = pp_paths
        self.output = output
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def write_input(self, atoms, properties, system_changes):
        """Write input parameters to files-file."""
        keys_to_silence = {"kpts", "kptopt", "ngkpt", "nshiftk", "shiftk"}
        # prevent to ASE to write kpoint section of the in File
        tmp_parameters = {p: self.parameters.get(p) for p in keys_to_silence if p in self.parameters}
        for param in tmp_parameters:
            if param == "kpts":
                self.parameters[param] = None
            del self.parameters[param]

        ase.io.abinit.write_all_inputs(
            atoms, properties, parameters=self.parameters,
            pp_paths=self.pp_paths,
            label=self.label, v8_legacy_format=self.v8_legacy_format)

        assert os.path.exists(f"{self.label}.in")
        # we have to fix the k-point input. ASE calculator can only create MP like meshes
        # we write this section ourselves
        buffer = io.StringIO()
        # restore parameters
        self.parameters.update(tmp_parameters)
        if tmp_parameters.get("kpts") is not None:
            kpoints = tmp_parameters.get("kpts")
            mp = kpts2mp(atoms, kpoints)
            buffer.write('kptopt %i\n' % (self.parameters.get("kptopt") or 1))
            buffer.write('ngkpt %d %d %d\n' % (self.parameters.get("ngkpt") or tuple(mp)))
            buffer.write('nshiftk %i\n' % (self.parameters.get("nshiftk") or 1))
            buffer.write('shiftk\n')
            buffer.write('%.3f %.3f %.3f\n' % (self.parameters.get("shiftk") or tuple((np.array(mp) + 1) % 2 * 0.5)))
            self.parameters["kpts"] = kpoints
            with open("%s.in" % self.label) as handle:
                buffer.write(handle.read())
            with open("%s.in" % self.label, 'w') as handle:
                handle.write(buffer.getvalue())

    def read(self, label):
        """Read results from ABINIT's text-output file."""
        # XXX I think we should redo the concept of 'restarting'.
        # It makes sense to load a previous calculation as
        #
        #  * static, calculator-independent results
        #  * an actual calculator capable of calculating
        #
        # Either of which is simpler than our current mechanism which
        # implies both at the same time.  Moreover, we don't need
        # something like calc.read(label).
        #
        # What we need for these two purposes is
        #
        #  * calc = MyCalculator.read(basefile)
        #      (or maybe it should return Atoms with calc attached)
        #  * results = read_results(basefile, format='abinit')
        #
        # where basefile determines the file tree.
        FileIOCalculator.read(self, label)
        self.atoms, self.parameters = ase.io.abinit.read_ase_and_abinit_inputs(self.label)
        self.results = ase.io.abinit.read_results(self.label, self._output_filename())

    def _output_filename(self):
        if self.v8_legacy_format:
            ext = '.txt'
        else:
            ext = '.abo'
        return self.label + ext

    def read_results(self):
        self.results = ase.io.abinit.read_results(self.label, self._output_filename())
        self.parameters.update(parse_abinit_input("%s.in" % self.label))

    def get_number_of_iterations(self):
        return self.results['niter']

    def get_electronic_temperature(self):
        return self.results['width']

    def get_number_of_electrons(self):
        return self.results['nelect']

    def get_number_of_bands(self):
        return self.results['nbands']

    def get_k_point_weights(self):
        return self.results['kpoint_weights']

    def get_bz_k_points(self):
        raise NotImplementedError

    def get_ibz_k_points(self):
        return self.results['ibz_kpoints']

    def get_spin_polarized(self):
        return self.results['eigenvalues'].shape[0] == 2

    def get_number_of_spins(self):
        return len(self.results['eigenvalues'])

    def get_fermi_level(self):
        return self.results['fermilevel']

    def get_eigenvalues(self, kpt=0, spin=0):
        return self.results['eigenvalues'][spin, kpt]

    def get_occupations(self, kpt=0, spin=0):
        raise NotImplementedError

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes, interval: float = 0.1):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise CalculatorSetupError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        command = self.command
        if 'PREFIX' in command:
            command = command.replace('PREFIX', self.prefix)

        if self.output == '-':
            stdout, stderr = (sys.stdout, sys.stderr)
            popen_kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            do_forward = True
            have_opened_file = False
        elif self.output is None:
            stdout, stderr = (os.devnull, os.devnull)
            popen_kwargs = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            do_forward = False
            have_opened_file = False
        elif isinstance(self.output, str):
            handle = open(f"{self.output}.out", 'w')
            stdout, stderr = (handle, handle)
            popen_kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            do_forward = True
            have_opened_file = True
        else:
            stdout, stderr = (self.output, self.output)
            popen_kwargs = dict(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            do_forward = True
            have_opened_file = False

        try:
            proc = subprocess.Popen(command, shell=True, cwd=self.directory, **popen_kwargs)
        except OSError as err:
            # Actually this may never happen with shell=True, since
            # probably the shell launches successfully.  But we soon want
            # to allow calling the subprocess directly, and then this
            # distinction (failed to launch vs failed to run) is useful.
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        if do_forward:
            stdout_src = io.TextIOWrapper(proc.stdout)
            stderr_src = io.TextIOWrapper(proc.stderr)
            while proc.poll() is None:
                stdout.write(stdout_src.read())
                stderr.write(stderr_src.read())
                time.sleep(interval)
            if have_opened_file:
                stdout.close()
                stderr.close()
        errorcode = proc.wait()

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                   '{} with error code {}'.format(self.name, command,
                                                  path, errorcode))
            raise CalculationFailed(msg)

        self.read_results()
