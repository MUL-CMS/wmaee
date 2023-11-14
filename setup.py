from setuptools import setup, find_packages

from _version import __version__

setup(
    name='wmaee',
    version=__version__,
    description='A collection of evaluation scripts for MUL-CMS group',
    author_email='david.holec@unileoben.ac.at',
    license='MIT',
    install_requires=['numpy', 'ase', 'matplotlib'],
    # install_requires=['numpy', 'scipy', 'ase', 'lammps', 'gpaw', 'kimpy', 'kim-query', 'matplotlib'],
    packages=find_packages('.', exclude=['tests', 'recipies']),
    include_package_data=True
)
