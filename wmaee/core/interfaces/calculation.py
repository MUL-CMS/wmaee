
import os
import abc
import tempfile
import contextlib
import dataclasses
from frozendict import frozendict
from typing import Any, TypeVar, Optional, Dict, Callable, NoReturn, TextIO, Union

from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.md.md import MolecularDynamics
from ase.calculators.calculator import Calculator
from ase.io.trajectory import TrajectoryReader, TrajectoryWriter

from wmaee.units import HARTREE_TO_EV
from wmaee.core.utils import merge, override_environ
from wmaee.core.interfaces.requirements import requires
from wmaee.core.interfaces.abinit_netcdf_trajectory import read_abinit_netcdf_trajectory
from wmaee.core.interfaces.runners import vasp, launch, write_input, read_results_vasp, gpaw, construct_calculator, abinit, read_results_abinit, calculate

Incar = Potcar = frozendict
MDCalculation = TypeVar("MDCalculation")
VaspCalculation = TypeVar("VaspCalculation")
GpawCalculation = TypeVar("GpawCalculation")
AbinitCalculation = TypeVar("AbinitCalculation")
GPAW = TypeVar("GPAW")


class Nothing:
    pass


def get_item(item: Any, *objects) -> Any:
    nothing = Nothing()
    for obj in objects:
        result = getattr(obj, item, nothing)
        if result is not nothing:
            return result
    raise AttributeError(item)


class AtomsAndCalculatorProxy(abc.ABC):

    def __getattr__(self, item):
        nothing = Nothing()
        result = self.__dict__.get(item, nothing)
        lookups = [self.atoms] + ([self.calculator] if self.calculator is not None else [])
        return get_item(item, self.atoms, *lookups) if result is nothing else result


class AbinitCalculation(AtomsAndCalculatorProxy):

    def __init__(self, atoms: Atoms, kpts=(1, 1, 1), ecut: float = 400, pps: str = "fhi", xc: str = "GGA", ixc: int = 11, prefix: Optional[str] = None, toldfe: float = 1e-4, **kwargs):
        self.atoms = atoms
        self.kpts = kpts
        self.pps = pps
        self.xc = xc
        self.prefix = self.atoms.get_chemical_formula() if prefix is None else prefix
        self.kwargs = frozendict(merge(dict(ecut=ecut, toldfe=toldfe, ixc=ixc, pps=pps, xc=xc), kwargs))
        self.calculator = None

    def set(self, *remove: str, **tags) -> AbinitCalculation:
        """
        set keyword arguments, override them or remove them

        :param remove: keyword argument to remove
        :type remove: str
        :param tags: keyword argument to set or override
        :return: calculation instance
        :rtype: VaspCalculation
        """
        new_kwargs = dict(self.kwargs, **tags)
        for removal in remove:
            if removal in new_kwargs:
                del new_kwargs[removal]
        self.kwargs = frozendict(new_kwargs)
        return self

    def run(self, directory: Optional[str] = None, ncpus: int = 2, output: Union[str, TextIO] = '-', **kwargs) -> NoReturn:
        """
        Run the AbinitCalculation instance. If {directory} is not `None` a temporary directory will be created.

        :param directory: the directory where the Abinit will be executed (default is `None`)
        :type directory: Optional[str]
        :param ncpus: number of MPI ranks (default is 2)
        :type ncpus: int
        :param output: redirection for the Abinit's output. "-" means stdout. Might be a filename of file object
            (default is "-")
        :type output: str | TextIO
        """
        directory = os.getcwd() if directory is None else directory
        output_kwargs = dict(ncpus=ncpus, directory=directory, output=output)
        abinit_kwargs = merge(output_kwargs, self.kwargs, kwargs)
        action = launch if self.calculator is None else calculate
        self.calculator, _ = abinit(action, self.atoms, prefix=self.prefix, kpts=self.kpts, **abinit_kwargs)

    @classmethod
    def from_directory(cls, path: Optional[str] = None, prefix: Optional[str] = None, **kwargs) -> AbinitCalculation:
        path = os.getcwd() if path is None else path
        if prefix is None:
            root, _, files = next(os.walk(path, topdown=True))
            files = set(files)
            input_file = next((f for f in files if f.endswith(".in")), None)
            if input_file is None:
                raise FileNotFoundError(f"No input file found in \"{root}\". Please specify a prefix!")
            prefix, *_ = input_file.split('.')

        calculator, atoms = read_results_abinit(prefix, directory=path, **kwargs)
        calculator.atoms = atoms

        calculation = cls(atoms, prefix=prefix, **kwargs)
        calculation.calculator = calculator
        return calculation

    def write_input(self, directory: Optional[str] = None, **kwargs) -> NoReturn:
        """
        Writes the input files the {directory}. {kwargs} are forwarded to the `ase.calculators.abinit.Abinit` constructor.

        :param directory: the directory to write the input files to (default is `os.getcwd()`)
        :type directory: str
        """
        directory = os.getcwd() if directory is None else directory
        if not os.path.exists(directory):
            os.makedirs(directory, exist_ok=True)
        abinit_kwargs = merge(self.kwargs, kwargs)
        self.calculator, _ = abinit(write_input, self.atoms, directory=directory, prefix=self.prefix, kpts=self.kpts, **abinit_kwargs)

    def get_trajectory(self):
        """
        Reads the file {label}o_HIST.nc and extracts the trajectory. Useful when visualizing e.g., relaxations
        """
        if self.calculator is None:
            raise RuntimeError("The calculation was not yet executed. Please call calculation.run() before loading the trajectory")
        return read_abinit_netcdf_trajectory(self.calculator)


@dataclasses.dataclass
class VaspCalculation(AtomsAndCalculatorProxy):
    """
    Interface to execute VASP calculations via `ase.calculators.vasp.Vasp`, with improved code semantics
    An example to use this class reads

        .. code-block:: python

            calculation = VaspCalculation(atoms, incar, kpts=(2,2,2), xc="pbe")
            calculation = calculation.tags(IBRION=2, LWAVE=False, LCHARG=True)  # override INCAR tags
            calculation.write_input("vasp_calc_dir_prepare")
            calculation.run()  # in temporary directory
            calculation.run("vasp_dir")  # run the calculation in the directory in "vasp_dir"

    """
    atoms: Atoms
    incar: frozendict = dataclasses.field(default_factory=Incar)
    kpts: Any = dataclasses.field(default=(1, 1, 1))
    potcar: frozendict = dataclasses.field(default_factory=lambda: frozendict(base="recommended"))
    gamma_centered: bool = dataclasses.field(default=True)
    xc: str = dataclasses.field(default="pbe")
    calculator: Optional[Vasp] = dataclasses.field(default=None)

    def tags(self, *remove: str, **tags) -> VaspCalculation:
        """
        set INCAR tags, override them or remove them

        :param remove: INCAR tags to remove
        :type remove: str
        :param tags: INCAR tags to set or override
        :return: calculation instance
        :rtype: VaspCalculation
        """
        new_incar = dict(self.incar, **tags)
        for removal in remove:
            if removal in new_incar:
                del new_incar[removal]
        self.incar = frozendict(new_incar)
        return self

    @classmethod
    def from_directory(cls, path: Optional[str] = None, **kwargs) -> VaspCalculation:
        """
        Create a VaspCalculation instance from a folder which contains VASP input files

        :param path: folder name to load the data from (default is `os.getcwd()`)
        :type path: str
        :param kwargs: keyword arguments are forwarded to `ase.calculator.vasp.Vasp` constructor
        :return: calculation instance
        :retype: VaspCalculation
        """
        path = os.getcwd() if path is None else path
        calculator, atoms = read_results_vasp(directory=path, **kwargs)
        # the VASP calculator class splits up the INCAR tags and gathers them in several different variables.
        # the names of these variables depend on the data type of the INCAR tags' value.
        parameters = ["int_params", "list_int_params", "float_params", "list_float_params", "bool_params",
                      "list_bool_params", "exp_params", "string_params", "special_params"]
        incar = frozendict({tag.upper(): value for parameter_collection in parameters for tag, value in
                            getattr(calculator, parameter_collection).items() if value is not None})
        return cls(atoms=calculator.get_atoms(), incar=incar, kpts=calculator.kpts,
                   gamma_centered=calculator.input_params.get("gamma", False), xc=calculator.get_xc_functional(),
                   calculator=calculator)

    def run(self, directory: Optional[str] = None, ncpus: int = 2, mode: str = "std", output: Union[str, TextIO] = '-', **kwargs) -> NoReturn:
        """
        Run the VaspCalculation instance. If {directory} is not None a temporary directory will be created.

        :param directory: the directory where the VASP code will be executed (default is `None`)
        :type directory: Optional[str]
        :param ncpus: number of MPI ranks (default is 2)
        :type ncpus: int
        :param mode: the VASP executable to execute. Either "std", "gam" or "ncl". (default is "std")
        :type mode: str
        :param output: redirection for the VASP's output. "-" means stdout. Might be a filename of file object
            (default is "-")
        :type output: str | TextIO
        """
        directory = os.getcwd() if directory is None else directory
        calc_kwargs = dict(kpts=self.kpts, xc=self.xc, setups=self.potcar, incar=self.incar,
                           gamma=self.gamma_centered)
        output_kwargs = dict(ncpus=ncpus, output=output, mode=mode, directory=directory)
        vasp_kwargs = merge(calc_kwargs, output_kwargs, **kwargs)
        action = launch if self.calculator is None else calculate
        self.calculator, _ = vasp(action, self.atoms, **vasp_kwargs)

    def write_input(self, directory: Optional[str] = None, **kwargs) -> NoReturn:
        """
        Writes the input files the {directory}. {kwargs} are forwarded to the `ase.calculators.vasp.Vasp` constructor.

        :param directory: the directory to write the input files to (default is `os.getcwd()`)
        :type directory: str
        """
        directory = os.getcwd() if directory is None else directory
        calc_kwargs = dict(kpts=self.kpts, xc=self.xc, setups=self.potcar, incar=self.incar,
                           gamma=self.gamma_centered)
        self.calculator, _ = vasp(write_input, self.atoms, **merge(calc_kwargs, dict(directory=directory), **kwargs))


class GpawCalculation(AtomsAndCalculatorProxy):

    def __init__(self, atoms: Atoms, kpts: Any = (1, 1, 1), xc: str = "PBE", **kwargs):
        """
        Initializes a GpawCalculation

        :param atoms: structure to compute
        :type atoms: ase.Atoms
        :param kpts: k-points object. {kpts} is passed on to the GPAW calculator constructor
        :param xc: the name of the XC functional to use (default is "PBE")
        :type xc: str
        :param kwargs: the keyword-arguments are forwarded to the GPAW calculator constructor
        """
        self.atoms = atoms
        self._input_parameters: Dict[str, Any] = merge(dict(kpts=kpts, xc=xc), **kwargs)
        self.calculator: Optional[GPAW] = None

    @requires("gpaw")
    def set(self, *remove:  str, **parameters) -> GpawCalculation:
        """
        set or remove parameter for a GPAW calculation. A parameter check is carried out with the new parameters, and
        an eventual exception if forwarded

        :param remove: parameters to remove
        :type remove: str
        :param parameters: name and value of the parameters for the GPAW calculation
        :return: calculation instance
        :rtype: GpawCalculation
        """
        tmp_parameters = merge(self._input_parameters, **parameters)
        for key in remove:
            if key in tmp_parameters:
                del tmp_parameters[key]
        try:
            from gpaw.new.ase_interface import InputParameters
            InputParameters(tmp_parameters)
        except:
            raise
        else:
            self._input_parameters = tmp_parameters
            return self

    @requires("gpaw")
    def run(self, directory: Optional[str] = None, ncpus: int = 2, output: Union[str, TextIO] = '-', **kwargs) -> NoReturn:
        """
        Executes the GPAW calculation.

        :param directory: the directory where the GPAW will be executed (default is `None`)
        :type directory: Optional[str]
        :param ncpus: number of threads (default is 2)
        :type ncpus: int
        :param output: redirection for the GPAW's output. "-" means stdout. Might be a filename of file object
            (default is "-")
        :type output: str | TextIO
        """
        directory_context = tempfile.TemporaryDirectory() if directory is None else contextlib.nullcontext(directory)
        with directory_context as wd:
            if self.calculator is None:
                output_kwargs = dict(ncpus=ncpus, output=output, directory=wd)
                gpaw_kwargs = merge(self._input_parameters, output_kwargs, **kwargs)
                self.calculator, _ = gpaw(launch, self.atoms, **gpaw_kwargs)
            else:
                launch(self.calculator, self.atoms)

    @requires("gpaw")
    def save(self, filename: str, mode: str = "") -> NoReturn:
        """
        Write the results of the calculation into a file. Traditionally it ends on .gpw.

        :param filename: the filename to write it to
        :type filename: str
        :param mode: if set to "all" wave functions will be included (default is "")
        :type mode: str
        """
        if self.calculator is None:
            raise RuntimeError("The calculation was not yet executed. Please call calculation.run() before saving it")
        self.calculator.write(filename, mode)

    @classmethod
    @requires("gpaw")
    def from_file(cls, filename: str) -> GpawCalculation:
         """
         Construct a GpawCalculation object from a file. Traditionally these files end with .gpw

         :param filename: the file to read from
         """
         calculator, _ = gpaw(construct_calculator, None, kpts=None, prefix="-", mode=None, xc=None, output=os.devnull)
         calculator.read(filename)
         calculation = GpawCalculation(calculator.atoms, kpts=None, xc=calculator.get_xc_functional())
         calculation.calculator = calculator
         calculation.atoms.calc = calculator
         return calculation


@dataclasses.dataclass
class MDCalculation(AtomsAndCalculatorProxy):
    atoms: Atoms
    model: str
    calculator: Optional[Calculator] = dataclasses.field(default=None)
    _dynamics: Optional[MolecularDynamics] = dataclasses.field(default=None)
    _filename: Optional[str] = dataclasses.field(default=None)
    _dump_interval: int = dataclasses.field(default=10)
    _trajectory_reader: Optional[TrajectoryReader] = dataclasses.field(default=None)

    def set(self, func, *args, **kwargs) -> MDCalculation:
        """
        Calls {func} and passes {atoms} as the first argument. Syntactic sugar for `func(self.atoms, *args, **kwargs)`.

        :param func: class of function to call
        :type func: Callable[[ase.Atoms, ...], NoReturn]
        :return: the calculation instance
        :rtype: MDCalculation
        """
        func(self.atoms, *args, **kwargs)
        return self

    @requires("kimpy")
    def dynamics(self, dynamics, *args, **kwargs) -> MDCalculation:
        """
        Calls {dynamics} and passes {atoms} as the first argument. Syntactic sugar for
        `dynamics(self.atoms, *args, **kwargs)`. Furthermore, it will set the internal dynamics object

        :param dynamics: class of function to call
        :type dynamics: Callable[[ase.Atoms, ...], NoReturn]
        :return: the calculation instance
        :rtype: MDCalculation
        """
        if self.calculator is None:
            from ase.calculators.kim import KIM
            self.calculator = KIM(self.model, options=dict(release_GIL=True))
            self.atoms.calc = self.calculator
        self._dynamics = dynamics(self.atoms, *args, **kwargs)
        return self

    def dump(self, filename: str, interval: int = 10) -> MDCalculation:
        if self._dynamics is None:
            raise RuntimeError("No dynamics was set yet. Use calculation.dynamics(...) before setting the dump")
        self._filename = filename
        self._dump_interval = interval
        return self

    def get_trajectory(self):
        if self._filename is None:
            raise RuntimeError("No trajectory filename was set")
        elif self._dynamics is None:
            raise RuntimeError("The calculation was not computed yet. Please calc calculation.run()")
        elif self._trajectory_reader is None:
            assert self._filename is not None
            self._trajectory_reader = TrajectoryReader(self._filename)
        return self._trajectory_reader

    def attach(self, f: Callable[[Atoms], NoReturn], interval: int = 50, pass_atoms: bool = False) -> MDCalculation:
        """
        Is a shortcut to `self._dynamics.attach(f)`. Is used to add callback functions. If {pass_atoms} is `True`
        {self.atoms} and {self._dynamics} will be passed to {f}

        :param f: class of function to call
        :type f: Callable[[ase.Atoms, ...], NoReturn]
        :param interval: number of steps after which {f} will be called
        :type interval: int
        :param pass_atoms: flag whether to pass the `ase.Atoms` and {dynamics} object to {f} (default is `False`)
        :type pass_atoms: bool
        :return: the calculation instance
        :rtype: MDCalculation
        """
        if self._dynamics is None:
            raise ValueError("No dynamics was defined for the MD calculation yet. "
                             "Use e.g. calculation.dynamics(Langevin, 5) to set the dynamics.")
        self._dynamics.attach((lambda: f(self.atoms, self._dynamics)) if pass_atoms else f, interval=interval)
        return self

    @requires("kimpy")
    def run(self, *args, ncpus: int = 2, **kwargs):
        if self._dynamics is None:
            raise ValueError("No dynamics was defined for the MD calculation yet. "
                             "Use e.g. calculation.dynamics(Langevin, 5) to set the dynamics.")
        with override_environ(OMP_NUM_THREADS=f"{ncpus}"):
            trajectory_context = TrajectoryWriter(self._filename, atoms=self.atoms) if self._filename is not None else contextlib.nullcontext(None)
            with trajectory_context as trajectory:
                if trajectory is not None:
                    self.attach(trajectory.write, interval=self._dump_interval)
                self._dynamics.run(*args, **kwargs)
        return self