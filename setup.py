from setuptools import setup, find_packages
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
