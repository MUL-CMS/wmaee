

from setuptools import setup, find_packages

setup(
    name='wmaee',
    version='0.1',
    description='A simple interface for the WMAAE exercises',
    author_email='dominik.gehringer@unileoben.ac.at',
    license='MIT',
    install_requires=['numpy', 'scipy', 'ase', 'lammps', 'gpaw', 'kimpy', 'kim-query', 'matplotlib'],
    packages=find_packages('.', exclude=['tests', 'recipies']),
    include_package_data=True
)