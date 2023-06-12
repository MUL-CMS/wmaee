import os
from ase import Atoms
from functools import partial
from frozendict import frozendict
from ase.calculators.calculator import Calculator
from typing import Optional, Any, Dict, NoReturn, Callable, Dict
from wmaee.core.utils import Config, render_command, override_environ

GAMMA_POINT = (1, 1, 1)

Command = Callable[[Calculator, Atoms], NoReturn]


def launch(_: Calculator, a: Atoms) -> NoReturn:
    a.get_total_energy()


def _vasp(action: Command, atoms: Atoms, kpts=GAMMA_POINT, xc="pbe", directory: str = os.getcwd(), ncpus: int = 2,
          version: str = "6.3.1", mode: str = "std", output: str = '-', setups: Dict = frozendict(base="recommended"),
          **kwargs):
    command = render_command("vasp", ncpus=ncpus, version=version, mode=mode)
    pseudo_potential_directory = Config().get("applications").get("vasp").get("potential_path")
    # override environment variables according to: https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html
    with override_environ(VASP_PP_PATH=pseudo_potential_directory):
        from ase.calculators.vasp import Vasp
        action(
            Vasp(atoms=atoms, directory=directory, command=command, txt=output, xc=xc, setups=setups, kpts=kpts,
                 **kwargs),
            atoms
        )
    return atoms


def _gpaw(action: Command, atoms: Atoms, kpts=GAMMA_POINT, mode: Optional[Any] = None, xc: str = "PBE",
                directory: str = os.getcwd(), ncpus: int = 2, prefix: Optional[str] = None, output: str = '-',
                **kwargs):
    application = "gpaw"
    prefix = prefix or atoms.get_chemical_formula()
    pseudo_potential_directory = Config().get("applications").get(application).get("potential_path")

    with override_environ(GPAW_SETUP_PATH=pseudo_potential_directory, OMP_NUM_THREADS=f"{ncpus}"):
        # It is crucial that the import is within the environmental variable override
        # the reason is that GPAW assembles the setup path not when constructing the GPAW but rather upon the module
        # import
        from gpaw import GPAW
        if mode is None:
            from gpaw import PW
            mode = PW()
        action(
            GPAW(atoms=atoms, mode=mode, xc=xc, directory=directory, txt=output, label=prefix, kpts=kpts, **kwargs),
            atoms
        )
    return atoms


launch_vasp = partial(_vasp, launch)
