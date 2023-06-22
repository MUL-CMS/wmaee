
"""
This code file is the additional "cost" of using Abinit. The ase.calculator.abinit implementation misses many features
as it has hardcoded parameters
"""

import io
import os
import re
import sys
import time
import warnings
import ase.units
import subprocess
import contextlib
import numpy as np
import ase.io.abinit
from frozendict import frozendict
from wmaee.core.utils import merge
from wmaee.units import HARTREE_TO_EV
from ase.calculators.abinit import Abinit
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


def _write_abinit_in_(fd, atoms, param=None, species=None, pseudos=None):
    import copy
    from ase.calculators.abinit import Abinit
    from ase.calculators.calculator import kpts2mp

    if param is None:
        param = {}

    _param = copy.deepcopy(Abinit.default_parameters)
    _param.update(param)
    param = _param

    if species is None:
        species = sorted(set(atoms.numbers))

    inp = {}
    inp.update(param)
    for key in ['xc', 'smearing', 'kpts', 'pps', 'raw']:
        del inp[key]

    smearing = param.get('smearing')
    if 'tsmear' in param or 'occopt' in param:
        assert smearing is None

    if smearing is not None:
        inp['occopt'] = {'fermi-dirac': 3,
                         'gaussian': 7}[smearing[0].lower()]
        inp['tsmear'] = smearing[1]

    inp['natom'] = len(atoms)

    if 'nbands' in param:
        inp['nband'] = param['nbands']
        del inp['nbands']

    # ixc is set from paw/xml file. Ignore 'xc' setting then.
    if param.get('pps') not in ['pawxml']:
        if 'ixc' not in param:
            inp['ixc'] = {'LDA': 7,
                          'PBE': 11,
                          'revPBE': 14,
                          'RPBE': 15,
                          'WC': 23}[param['xc']]

    magmoms = atoms.get_initial_magnetic_moments()
    if magmoms.any():
        inp['nsppol'] = 2
        fd.write('spinat\n')
        for n, M in enumerate(magmoms):
            fd.write('%.14f %.14f %.14f\n' % (0, 0, M))
    else:
        inp['nsppol'] = 1

    if param['kpts'] is not None:
        mp = kpts2mp(atoms, param['kpts'])
        fd.write('kptopt %i\n' % (param.get("kptopt") or 1))
        fd.write('ngkpt %d %d %d\n' % (param.get("ngkpt") or tuple(mp)))
        fd.write('nshiftk %i\n' % (param.get("nshiftk") or 1))
        fd.write('shiftk\n')
        fd.write('%.3f %.3f %.3f\n' % (param.get("shiftk") or tuple((np.array(mp) + 1) % 2 * 0.5)))

    valid_lists = (list, np.ndarray)
    for key in sorted(inp):
        value = inp[key]
        unit = ase.io.abinit.keys_with_units.get(key)
        if unit is not None:
            if 'fs**2' in unit:
                value /= ase.units.fs ** 2
            elif 'fs' in unit:
                value /= ase.units.fs
        if isinstance(value, valid_lists):
            if isinstance(value[0], valid_lists):
                fd.write("{}\n".format(key))
                for dim in value:
                    ase.io.abinit.write_list(fd, dim, unit)
            else:
                fd.write("{}\n".format(key))
                ase.io.abinit.write_list(fd, value, unit)
        else:
            if unit is None:
                fd.write("{} {}\n".format(key, value))
            else:
                fd.write("{} {} {}\n".format(key, value, unit))

        if param['raw'] is not None:
            if isinstance(param['raw'], str):
                raise TypeError('The raw parameter is a single string; expected '
                                'a sequence of lines')
            for line in param['raw']:
                if isinstance(line, tuple):
                    fd.write(' '.join(['%s' % x for x in line]) + '\n')
                else:
                    fd.write('%s\n' % line)

        fd.write('#Definition of the unit cell\n')
        fd.write('acell\n')
        fd.write('%.14f %.14f %.14f Angstrom\n' % (1.0, 1.0, 1.0))
        fd.write('rprim\n')
        if atoms.cell.rank != 3:
            raise RuntimeError('Abinit requires a 3D cell, but cell is {}'
                               .format(atoms.cell))
        for v in atoms.cell:
            fd.write('%.14f %.14f %.14f\n' % tuple(v))

        fd.write('chkprim %i # Allow non-primitive cells\n' % (param.get("chkprim") or 0))

        fd.write('#Definition of the atom types\n')
        fd.write('ntypat %d\n' % (len(species)))
        fd.write('znucl {}\n'.format(' '.join(str(Z) for Z in species)))
        fd.write('#Enumerate different atomic species\n')
        fd.write('typat')
        fd.write('\n')

        types = []
        for Z in atoms.numbers:
            for n, Zs in enumerate(species):
                if Z == Zs:
                    types.append(n + 1)
        n_entries_int = 20  # integer entries per line
        for n, type in enumerate(types):
            fd.write(' %d' % (type))
            if n > 1 and ((n % n_entries_int) == 1):
                fd.write('\n')
        fd.write('\n')

        if pseudos is not None:
            listing = ',\n'.join(pseudos)
            line = f'pseudos "{listing}"\n'
            fd.write(line)

        fd.write('#Definition of the atoms\n')
        fd.write('xcart\n')
        for pos in atoms.positions / ase.units.Bohr:
            fd.write('%.14f %.14f %.14f\n' % tuple(pos))

        fd.write('chkexit 1 # abinit.exit file in the running '
                 'directory terminates after the current SCF\n')


def write_abinit_in(fd, atoms, param=None, species=None, pseudos=None):
    import copy
    from ase.calculators.calculator import kpts2mp
    from ase.calculators.abinit import Abinit

    written_properties = set()
    if param is None:
        param = {}

    _param = copy.deepcopy(Abinit.default_parameters)
    _param.update(param)
    param = _param

    if species is None:
        species = sorted(set(atoms.numbers))

    inp = {}
    inp.update(param)
    for key in ['xc', 'smearing', 'kpts', 'pps', 'raw']:
        del inp[key]

    smearing = param.get('smearing')
    if 'tsmear' in param or 'occopt' in param:
        assert smearing is None

    if smearing is not None:
        inp['occopt'] = {'fermi-dirac': 3,
                         'gaussian': 7}[smearing[0].lower()]
        inp['tsmear'] = smearing[1]

    inp['natom'] = len(atoms)

    if 'nbands' in param:
        inp['nband'] = param['nbands']
        del inp['nbands']

    # ixc is set from paw/xml file. Ignore 'xc' setting then.
    if param.get('pps') not in ['pawxml']:
        if 'ixc' not in param:
            inp['ixc'] = {'LDA': 7,
                          'PBE': 11,
                          'revPBE': 14,
                          'RPBE': 15,
                          'WC': 23}[param['xc']]

    magmoms = atoms.get_initial_magnetic_moments()
    if magmoms.any():
        inp['nsppol'] = 2
        fd.write('spinat\n')
        for n, M in enumerate(magmoms):
            fd.write('%.14f %.14f %.14f\n' % (0, 0, M))
        written_properties |= {"nsppol", "spinat"}
    else:
        inp['nsppol'] = 1
        written_properties.add("nsppol")

    if param['kpts'] is not None:
        mp = kpts2mp(atoms, param['kpts'])
        fd.write('kptopt %i\n' % (param.get("kptopt") or 1))
        fd.write('ngkpt %d %d %d\n' % (param.get("ngkpt") or tuple(mp)))
        fd.write('nshiftk %i\n' % (param.get("nshiftk") or 1))
        fd.write('shiftk\n')
        fd.write('%.3f %.3f %.3f\n' % (param.get("shiftk") or tuple((np.array(mp) + 1) % 2 * 0.5)))
        written_properties |= {"kptopt", "ngkpt", "nshiftk", "shiftk"}

    fd.write('chkprim %i # Allow non-primitive cells\n' % (param.get("chkprim") or 0))
    written_properties.add("chkprim")

    valid_lists = (list, np.ndarray)
    for key in sorted(inp):
        if key in written_properties:
            continue
        value = inp[key]
        unit = ase.io.abinit.keys_with_units.get(key)
        if unit is not None:
            if 'fs**2' in unit:
                value /= ase.units.fs**2
            elif 'fs' in unit:
                value /= ase.units.fs

        if isinstance(value, valid_lists):
            if isinstance(value[0], valid_lists):
                fd.write("{}\n".format(key))
                for dim in value:
                    ase.io.abinit.write_list(fd, dim, unit)
            else:
                fd.write("{}\n".format(key))
                ase.io.abinit.write_list(fd, value, unit)
        else:
            if unit is None:
                fd.write("{} {}\n".format(key, value))
            else:
                fd.write("{} {} {}\n".format(key, value, unit))

    if param['raw'] is not None:
        if isinstance(param['raw'], str):
            raise TypeError('The raw parameter is a single string; expected '
                            'a sequence of lines')
        for line in param['raw']:
            if isinstance(line, tuple):
                fd.write(' '.join(['%s' % x for x in line]) + '\n')
            else:
                fd.write('%s\n' % line)

    fd.write('#Definition of the unit cell\n')
    fd.write('acell\n')
    fd.write('%.14f %.14f %.14f Angstrom\n' % (1.0, 1.0, 1.0))
    fd.write('rprim\n')
    if atoms.cell.rank != 3:
        raise RuntimeError('Abinit requires a 3D cell, but cell is {}'
                           .format(atoms.cell))
    for v in atoms.cell:
        fd.write('%.14f %.14f %.14f\n' % tuple(v))


    fd.write('#Definition of the atom types\n')
    fd.write('ntypat %d\n' % (len(species)))
    fd.write('znucl {}\n'.format(' '.join(str(Z) for Z in species)))
    fd.write('#Enumerate different atomic species\n')
    fd.write('typat')
    fd.write('\n')

    types = []
    for Z in atoms.numbers:
        for n, Zs in enumerate(species):
            if Z == Zs:
                types.append(n + 1)
    n_entries_int = 20  # integer entries per line
    for n, type in enumerate(types):
        fd.write(' %d' % (type))
        if n > 1 and ((n % n_entries_int) == 1):
            fd.write('\n')
    fd.write('\n')

    if pseudos is not None:
        listing = ',\n'.join(pseudos)
        line = f'pseudos "{listing}"\n'
        fd.write(line)

    fd.write('#Definition of the atoms\n')
    fd.write('xcart\n')
    for pos in atoms.positions / ase.units.Bohr:
        fd.write('%.14f %.14f %.14f\n' % tuple(pos))

    fd.write('chkexit 1 # abinit.exit file in the running '
             'directory terminates after the current SCF\n')


def write_all_inputs(atoms, properties, parameters,
                     pp_paths=None,
                     raise_exception=True,
                     label='abinit',
                     *, v8_legacy_format=True):
    species = sorted(set(atoms.numbers))
    if pp_paths is None:
        pp_paths = ase.io.abinit.get_default_abinit_pp_paths()
    ppp = ase.io.abinit.get_ppp_list(atoms, species,
                       raise_exception=raise_exception,
                       xc=parameters.xc,
                       pps=parameters.pps,
                       search_paths=pp_paths)

    if v8_legacy_format is None:
        warnings.warn(ase.io.abinit.abinit_input_version_warning,
                      FutureWarning)
        v8_legacy_format = True

    if v8_legacy_format:
        with open(label + '.files', 'w') as fd:
            ase.io.abinit.write_files_file(fd, label, ppp)
        pseudos = None

        # XXX here we build the txt filename again, which is bad
        # (also defined in the calculator)
        output_filename = label + '.txt'
    else:
        pseudos = ppp  # Include pseudopotentials in inputfile
        output_filename = label + '.abo'

    # Abinit will write to label.txtA if label.txt already exists,
    # so we remove it if it's there:
    if os.path.isfile(output_filename):
        os.remove(output_filename)

    parameters.write(label + '.ase')

    with open(label + '.in', 'w') as fd:
        write_abinit_in(fd, atoms, param=parameters, species=species,
                        pseudos=pseudos)


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
        write_all_inputs(
            atoms, properties, parameters=self.parameters,
            pp_paths=self.pp_paths,
            label=self.label, v8_legacy_format=self.v8_legacy_format)

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
