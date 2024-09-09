import glob
import multiprocessing
import os
import abc
import shutil
import unittest
import functools
import numpy as np
from typing import Optional
import matplotlib.pyplot as plt
from numpy.typing import ArrayLike
from wmaee.core.utils import override_environ
from wmaee import Atoms, GpawCalculation, AbinitCalculation
from wmaee.extensions.energy_volume import birch_murnaghan_fit


AL_KPTS = (4, 4, 4)
AL_EQ = 4.05
AL_MIN = 0.95
AL_MAX = 1.05


def aluminium(a: float) -> Atoms:
    return Atoms(cell=np.eye(3) * a, symbols=['Al'] * 4,
                 scaled_positions=[[0.5, 0.5, 0], [0.5, 0, 0.5, ], [0, 0.5, 0.5], [0.0, 0.0, 0.0]],
                 pbc=True)


def monkey_patch():
    import ase.utils
    if not hasattr(ase.utils, "StringIO"):
        import io
        ase.utils.StringIO = io.StringIO


def inject_config_file(f):

    @functools.wraps(f)
    def _wrapper(*args, **kwargs):
        this_folder = os.path.dirname(__file__)
        config_file = os.path.join(this_folder, "wmaee.conf.yaml")
        with override_environ(WMAEE_CONFIG_FILE=config_file):
            return f(*args, **kwargs)

    return _wrapper


def fit_and_plot(volumes: ArrayLike, energies: ArrayLike, filename: str, title: Optional[str] = None):
    (e0, v0, b0, bp), fit = birch_murnaghan_fit(volumes, energies)
    vmin, vmax = np.amin(volumes), np.amax(volumes)
    vaxis = np.linspace(vmin, vmax)
    plt.plot(vaxis, fit(vaxis), color='r', label="Birch Murnaghan EOS - fit", zorder=-10)
    if title is not None:
        plt.title(title)
    plt.scatter(volumes, energies)
    plt.axvline(v0, color='b', linestyle='--')
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.text(0.75*(xmax - xmin) + xmin, 0.9*(ymax - ymin) + ymin, f"$a_0 = {v0**(1/3):.2f} A$", verticalalignment='center', horizontalalignment='center')
    plt.text(0.75*(xmax - xmin) + xmin, 0.85*(ymax - ymin) + ymin, f"$V_0 = {v0:.2f} A^3$", verticalalignment='center', horizontalalignment='center')
    plt.text(0.75*(xmax - xmin) + xmin, 0.8*(ymax - ymin) + ymin, f"$E_0 = {e0:.2f}$ eV", verticalalignment='center', horizontalalignment='center')
    plt.text(0.75*(xmax - xmin) + xmin, 0.75*(ymax - ymin) + ymin, f"$B_0 = {b0*160.19:.2f}$ GPa", verticalalignment='center', horizontalalignment='center')
    plt.savefig(filename)
    return e0, v0, b0, bp


class TestCalculation:
    def setup(self) -> None:
        monkey_patch()
        self.a_eq = AL_EQ
        self.a_min = AL_MIN
        self.a_max = AL_MAX
        self.ncalc = 6
        self.kpts = AL_KPTS
        self.this_folder = os.path.dirname(__file__)

    @inject_config_file
    def compute_eos(self, plotfile, title: Optional[str] = None):
        volumes, energies = [], []
        for a in np.linspace(self.a_min, self.a_max, self.ncalc):
            filename = self.calculation_filename(a)
            if not os.path.exists(filename):
                atoms = aluminium(a * self.a_eq)
                calculation = self.make_calculation(atoms, self.kpts)
                self.calculate(calculation, filename)
                self.save(calculation, filename)
            else:
                calculation = self.from_disk(filename)
                atoms = calculation.atoms
            volumes.append(atoms.cell.volume)
            energies.append(atoms.get_total_energy())
        return fit_and_plot(volumes, energies, os.path.join(self.this_folder, plotfile), title=title)

    @inject_config_file
    def load_eos(self, plotfile, title: Optional[str] = None):
        volumes, energies = [], []
        for a in np.linspace(self.a_min, self.a_max, self.ncalc):
            filename = self.calculation_filename(a)
            calculation = self.from_disk(filename)
            atoms = calculation.atoms
            volumes.append(atoms.cell.volume)
            energies.append(atoms.get_total_energy())
        return fit_and_plot(volumes, energies, os.path.join(self.this_folder, plotfile), title=title)

    @abc.abstractmethod
    def calculation_filename(self, a: float):
        raise NotImplementedError

    @abc.abstractmethod
    def make_calculation(self, atoms, kpts):
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, calculation, filename):
        raise NotImplementedError

    @abc.abstractmethod
    def from_disk(self, path):
        raise NotImplementedError

    @abc.abstractmethod
    def calculate(self, calculation, filename):
        raise NotImplementedError


class TestGpawCalculation(unittest.TestCase, TestCalculation):
    def setUp(self) -> None:
        self.setup()

    def calculation_filename(self, a: float):
        return os.path.join(self.this_folder, f"a-{a:.2f}.gpw")

    def make_calculation(self, atoms, kpts):
        return GpawCalculation(atoms, kpts=kpts)

    def save(self, calculation, filename):
        calculation.save(filename)

    def calculate(self, calculation, _):
        calculation.run()

    def from_disk(self, path):
        return GpawCalculation.from_file(path)

    @inject_config_file
    def test_eos(self):
        e0, v0, b0, bp = self.compute_eos("eos-bm-gpaw-calc.pdf", title="GPAW - $k$-points: $\Gamma$ - (4,4,4)")
        self.assertTrue(4.04 < v0 ** (1/3) < 4.06)

    @inject_config_file
    def test_load_eos(self):
        e0, v0, b0, bp = self.load_eos("eos-bm-gpaw-load.pdf", title="GPAW - $k$-points: $\Gamma$ - (4,4,4)")
        self.assertTrue(4.04 < v0 ** (1 / 3) < 4.06)

    def tearDownClass(*_) -> None:
        for gpwfile in glob.glob(os.path.join(os.path.dirname(__file__), "*.gpw")):
            os.remove(gpwfile)


class TestAbinitCalculation(unittest.TestCase, TestCalculation):
    def setUp(self) -> None:
        self.setup()

    def calculation_filename(self, a: float):
        return os.path.join(self.this_folder, f"abinit-a-{a:.2f}")

    def make_calculation(self, atoms, kpts):
        return AbinitCalculation(atoms, kpts=kpts)

    def save(self, calculation, filename):
        pass

    def calculate(self, calculation, directory):
        calculation.run(directory=directory, ncpus=multiprocessing.cpu_count())

    def from_disk(self, path):
        return AbinitCalculation.from_directory(path)

    @inject_config_file
    def test_eos(self):
        e0, v0, b0, bp = self.compute_eos("eos-bm-abinit-calc.pdf", title="Abinit - $k$-points: $\Gamma$ - (4,4,4)")
        self.assertTrue(4.07 < v0 ** (1/3) < 4.09)

    @inject_config_file
    def test_load_eos(self):
        e0, v0, b0, bp = self.load_eos("eos-bm-abinit-load.pdf", title="Abinit - $k$-points: $\Gamma$ - (4,4,4)")
        self.assertTrue(4.07 < v0 ** (1 / 3) < 4.09)

    def tearDownClass(*_) -> None:
        for abinit_dir in glob.glob(os.path.join(os.path.dirname(__file__), "abinit-a*")):
            shutil.rmtree(abinit_dir)

if __name__ == "__main__":
    unittest.main()