from setuptools import setup, find_packages

setup(
    name='wmaee',
    version='2.2',
    description='A collection of evaluation scripts for MUL-CMS group',
    author_email='david.holec@unileoben.ac.at',
    license='MIT',
    install_requires=['numpy', 'ase', 'matplotlib'],
    # install_requires=['numpy', 'scipy', 'ase', 'lammps', 'gpaw', 'kimpy', 'kim-query', 'matplotlib'],
    packages=find_packages('.', exclude=['tests', 'recipies']),
    include_package_data=True
)
