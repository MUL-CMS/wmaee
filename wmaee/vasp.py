import logging
from wmaee.core.common import working_directory, LoggerMixin
from wmaee.utils import ase_to_pymatgen
from wmaee.potentials import construct_potcar
from wmaee.core.runner import vasp
from pymatgen.io.vasp import Poscar, Potcar, Kpoints, Incar, Vasprun
from pymatgen import Structure
from ase import Atoms
from uuid import uuid4
from os import getcwd
from os.path import exists

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
        _write_input(self._structure, self._incar, kpoints=self._kpoints, potcar=self._potcar, xc_func=self._xc_func,
                     directory=directory)

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
    _write_input(inp.structure, inp.incar, kpoints=inp.kpoints, potcar=inp.potcar, xc_func=inp.xc_func,
                 directory=directory)


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