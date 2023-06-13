
import os
import tempfile
import contextlib
import dataclasses
from frozendict import frozendict
from typing import NoReturn, Any, TypeVar, Optional, Dict

from ase import Atoms
from ase.calculators.vasp import Vasp

from wmaee.core.utils import merge
from wmaee.core.interfaces.requirements import requires
from wmaee.core.interfaces.runners import vasp, launch, write_input, read_results_vasp, gpaw, construct_calculator

Incar = Potcar = frozendict
VaspCalculation = TypeVar("VaspCalculation")
GPAW = TypeVar("GPAW")


@dataclasses.dataclass
class VaspCalculation:
    atoms: Atoms
    incar: frozendict = dataclasses.field(default_factory=Incar)
    kpts: Any = dataclasses.field(default=(1, 1, 1))
    potcar: frozendict = dataclasses.field(default=frozendict(base="recommended"))
    gamma_centered: bool = dataclasses.field(default=True)
    xc: str = dataclasses.field(default="pbe")
    calculator: Optional[Vasp] = dataclasses.field(default=None)

    def tags(self, *remove: str, **tags) -> VaspCalculation:
        new_incar = dict(self.incar, **tags)
        for removal in remove:
            if removal in new_incar:
                del new_incar[removal]
        self.incar = frozendict(new_incar)
        return self

    @classmethod
    def from_directory(cls, path: str = os.getcwd(), **kwargs) -> VaspCalculation:
        calculator, atoms = read_results_vasp(directory=path, **kwargs)
        # the vasp calculator class splits up the incar tags and gathers them in several different variables.
        # the names of these variables depend on the data type of the incar tags' value.
        parameters = ["int_params", "list_int_params", "float_params", "list_float_params", "bool_params",
                      "list_bool_params", "exp_params", "string_params", "special_params"]
        incar = frozendict({tag.upper(): value for parameter_collection in parameters for tag, value in
                            getattr(calculator, parameter_collection).items() if value is not None})
        return cls(atoms=calculator.get_atoms(), incar=incar, kpts=calculator.kpts,
                   gamma_centered=calculator.input_params.get("gamma", False), xc=calculator.get_xc_functional(),
                   calculator=calculator)

    def run(self, directory: Optional[str] = None, ncpus: int = 2, mode: str = "std", output: str = '-', **kwargs):
        directory_context = tempfile.TemporaryDirectory() if directory is None else contextlib.nullcontext(directory)

        with directory_context as wd:
            calc_kwargs = dict(kpts=self.kpts, xc=self.xc, setups=self.potcar, incar=self.incar,
                               gamma=self.gamma_centered)
            output_kwargs = dict(ncpus=ncpus, output=output, mode=mode, directory=wd)
            vasp_kwargs = merge(calc_kwargs, output_kwargs, **kwargs)
            self.calculator, _ = vasp(launch, self.atoms, **vasp_kwargs)

    def write_input(self, directory: str = os.getcwd(), **kwargs):
        calc_kwargs = dict(kpts=self.kpts, xc=self.xc, setups=self.potcar, incar=self.incar,
                           gamma=self.gamma_centered)
        _ = vasp(write_input, self.atoms, **merge(calc_kwargs, dict(directory=directory), **kwargs))


class GpawCalculation:

    def __init__(self, atoms: Atoms, kpts: Any = (1, 1, 1), xc: str = "PBE", **kwargs):
        self.atoms = atoms
        self._input_parameters: Dict[str, Any] = merge(dict(kpts=kpts, xc=xc), **kwargs)
        self.calculator: Optional[GPAW] = None

    @requires("gpaw")
    def set(self, *remove, **kwargs):
        tmp_parameters = merge(self._input_parameters, **kwargs)
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
    def run(self, directory: Optional[str] = None, ncpus: int = 2, output: str = '-', **kwargs):
        directory_context = tempfile.TemporaryDirectory() if directory is None else contextlib.nullcontext(directory)

        with directory_context as wd:
            output_kwargs = dict(ncpus=ncpus, output=output, directory=wd)
            gpaw_kwargs = merge(self._input_parameters, output_kwargs, **kwargs)
            self.calculator, _ = gpaw(launch, self.atoms, **gpaw_kwargs)


@dataclasses.dataclass
class MDCalculation:
    atoms: Atoms
    model: str
    _calculator: Optional[Calculator] = dataclasses.field(default=None)
    _dynamics: Optional[MolecularDynamics] = dataclasses.field(default=None)

    def set(self, f, *args, **kwargs):
        f(self.atoms, *args, **kwargs)
        return self

    def dynamics(self, dynamics, *args, **kwargs):
        if self._calculator is None:
            from ase.calculators.kim import KIM
            self._calculator = KIM(self.model)
            self.atoms.calc = self._calculator
        self._dynamics = dynamics(self.atoms, *args, **kwargs)
        return self

    def attach(self, f: Callable[[Atoms, ParamSpecKwargs], NoReturn], interval: int = 50, pass_atoms: bool = False):
        if self._dynamics is None:
            raise ValueError("No dynamics was defined for the MD calculation yet. "
                             "Use e.g. calculation.dynamics(Langevin, 5) to set the dynamics.")
        self._dynamics.attach((lambda: f(self.atoms, self._dynamics)) if pass_atoms else f, interval=interval)
        return self

    @requires("kimpy")
    def run(self, *args, **kwargs):
        if self._dynamics is None:
            raise ValueError("No dynamics was defined for the MD calculation yet. "
                             "Use e.g. calculation.dynamics(Langevin, 5) to set the dynamics.")
        self._dynamics.run(*args, **kwargs)
        return self