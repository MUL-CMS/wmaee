from setuptools import setup, find_packages
from setuptools_rust import RustExtension, Binding
from pkg_resources import resource_filename
from fnmatch import fnmatch
from os.path import isfile, join
from os import listdir

__author__ = "Dominik Gehringer"
__email__ = "dgehringer@protonmail.com"
__maintainer__ = "Dominik Gehringer"
__license__ = "New BSD License"


setup(
    name="wmaee",
    version="0.1",
    author='Dominik Gehringer',
    author_email='dominik.gehringer@unileoben.ac.at',
    packages=find_packages(),
    # package_data=package_data,
    rust_extensions = [
        RustExtension("wmaee.extensions.external.parse_vasp_volumetric",
                      binding=Binding.PyO3,
                      path="wmaee/extensions/external/parse-vasp-volumetric-rs/Cargo.toml")
    ],
    include_package_data=True,
    install_requires=[
        "pymatgen",
        "ase",
        "plotly",
        "matplotlib",
        "numpy",
        "scipy",
        "pyiron",
        "nglview"
    ],
    python_requires='>=3.7'
)
