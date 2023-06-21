
import os
import sys
import tempfile
import contextlib
from frozendict import frozendict
from typing import Optional, Any, Callable, Dict, Tuple, Iterable

from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.calculators.abinit import Abinit
from ase.calculators.calculator import Calculator, all_changes

from wmaee.core.interfaces.requirements import requires
from wmaee.core.interfaces.utils import Config, render_command
from wmaee.core.utils import override_environ, working_directory

GAMMA_POINT = (1, 1, 1)



ABI_UNIT_NAMES = {
    s.lower() for s in (
        "au", "nm",
        "Angstr", "Angstrom", "Angstroms", "Bohr", "Bohrs",
        "eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs",
        "T", "Tesla",)
}

Command = Callable[[Calculator, Atoms], Any]


def launch(c: Calculator, a: Atoms) -> Any:
    a.get_total_energy()
    return c, a


def construct_calculator(c: Calculator, a: Atoms) -> Any:
    return c, a


def write_input(c: Calculator, a: Atoms) -> Any:
    c.write_input(a, ("energies",), tuple(all_changes))
    return c, a


def calculate(c: Calculator, a: Atoms) -> Any:
    c.calculate(a, ("energies",), tuple(all_changes))
    return c, a


def vasp(action: Command, atoms: Atoms, kpts=GAMMA_POINT, xc="pbe", directory: Optional[str] = None, ncpus: int = 2,
         version: str = "6.3.1", mode: str = "std", output: str = '-', setups: Dict = frozendict(base="recommended"),
         incar: Dict[str, Any] = frozendict(), gamma: bool = True, **kwargs):
    command = render_command("vasp", ncpus=ncpus, version=version, mode=mode)
    pseudo_potential_directory = Config().get("applications").get("vasp").get("potential_path")
    # override environment variables according to: https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html
    directory = os.getcwd() if directory is None else directory
    with override_environ(VASP_PP_PATH=pseudo_potential_directory):
        from ase.calculators.vasp import Vasp
        calculator = Vasp(
            atoms=atoms, directory=directory, command=command, txt=output, xc=xc, setups=setups, kpts=kpts,
            gamma=gamma, **kwargs
        )
        calculator.set(**{incar_key.lower(): value for incar_key, value in incar.items()})
        return action(calculator, atoms)


def abinit(action: Command, atoms: Atoms, kpts=GAMMA_POINT, pps: str = "fhi", xc: str = "LDA", ncpus: int = 2, prefix: Optional[str] = None,  output: str = '-', directory: Optional[str] = None, **kwargs):
    directory = os.getcwd() if directory is None else directory
    prefix = atoms.get_chemical_formula() if prefix is None else prefix
    command = render_command("abinit", ncpus=ncpus, prefix=prefix)
    pseudo_potential_directory = Config().get("applications").get("abinit").get("potential_path")
    subdirectories = frozenset({"LDA_FHI", "GGA_FHI", "LDA_HGH", "LDA_PAW", "LDA_TM", "GGA_FHI", "GGA_HGHK", "GGA_PAW"})
    pseudo_potential_path = ":".join(os.path.join(pseudo_potential_directory, subdir) for subdir in subdirectories)

    if output == '-':
        output_file = sys.stdout
    elif isinstance(output, str):
        output_file = open(output)
    else:
        output_file = output

    with override_environ(ABINIT_PP_PATH=pseudo_potential_path, ASE_ABINIT_COMMAND=command):
        from ase.calculators.abinit import Abinit
        with contextlib.redirect_stdout(output_file):
            return action(
                Abinit(atoms=atoms, label=prefix, pps=pps, xc=xc, directory=directory, kpts=kpts, v8_legacy_format=False, **kwargs),
                atoms
            )


@requires("gpaw")
def gpaw(action: Command, atoms: Atoms, kpts=GAMMA_POINT, mode: Optional[Any] = None, xc: str = "PBE",
         directory: str = os.getcwd(), ncpus: int = 2, prefix: Optional[str] = None, output: str = '-',
         **kwargs):
    application = "gpaw"
    prefix = prefix or atoms.get_chemical_formula()
    pseudo_potential_directory = Config().get("applications").get(application).get("potential_path")

    with override_environ(GPAW_SETUP_PATH=pseudo_potential_directory, OMP_NUM_THREADS=f"{ncpus}"):

        # It is crucial that the import is within the environmental variable override
        # the reason is that GPAW assembles the setup path not when constructing the GPAW but rather upon the module
        # import
        import gpaw
        from gpaw import GPAW
        if pseudo_potential_directory not in gpaw.setup_paths:
            gpaw.setup_paths.append(pseudo_potential_directory)
        if mode is None:
            from gpaw import PW
            mode = PW()
        return action(
            GPAW(atoms=atoms, mode=mode, xc=xc, directory=directory, txt=output, label=prefix, kpts=kpts, **kwargs),
            atoms
        )


def read_results_vasp(directory: Optional[str] = None, **kwargs) -> Tuple[Vasp, Atoms]:
    """
    Create a `ase.calculator.vasp.Vasp` instance, by reading the contents of {directory}

    :param directory: the directory to read from
    :type directory: str
    :param kwargs: are forwarded to `ase.calculator.vasp.Vasp` constructor
    :return: the `ase.calculator.vasp.Vasp` calculator instance and the `ase.Atom` object
    :rtype: Tuple[Vasp, Atoms]:
    """
    directory = os.getcwd() if directory is None else directory
    calculator = Vasp(directory=directory, restart=True, **kwargs)
    calculator.read_results()
    return calculator, calculator.get_atoms()


def read_results_abinit(prefix: str, directory: Optional[str] = None, **kwargs) -> Tuple[Abinit, Atoms]:
    """
    Create a `ase.calculators.abinit.Abinit` instance, by reading the contents of {directory}

    :param prefix: abinit prefix for the Abinit input and output files
    :type prefix: str
    :param directory: the directory to read from
    :type directory: str
    :param kwargs: are forwarded to `ase.calculators.abinit.Abinit` constructor
    :return: the `ase.calculators.abinit.Abinit` calculator instance and the `ase.Atom` object
    :rtype: Tuple[Abinit, Atoms]:
    """

    directory = os.getcwd() if directory is None else directory
    with working_directory(directory):
        calculator = Abinit(directory=directory, restart=prefix, label=prefix, v8_legacy_format=False, **kwargs)
    calculator.read_results()
    return calculator, calculator.get_atoms()

