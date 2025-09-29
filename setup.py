from setuptools import setup, find_packages

setup(
    name='wmaee',
    version='2.3',
    description='A collection of evaluation scripts for MUL-CMS group',
    author='David Holec',
    author_email='david.holec@unileoben.ac.at',
    license='MIT',
    packages=find_packages('.', exclude=['tests', 'recipies']),
    include_package_data=True,
    install_requires=[
        'numpy',
        'ase',
        'matplotlib',
        'scipy'
    ],
    extras_require={
        'pyiron': [
            'pyiron_atomistics',
            'pyiron_base'  
        ],
        'atomistics': [
            'lammps',
            'kimpy',
            'kim-query'
        ],
        'dft': [
            'gpaw'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    long_description=(
        'Some parts of this package may require external software '
        'such as VASP or LAMMPS for full functionality.'
    ),
)
