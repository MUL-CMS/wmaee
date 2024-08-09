"""
==============================================================================
Module Name: pyiron_NEB_task.py
==============================================================================
Description:
    This module provides a framework for setting up and executing NEB (Nudged
    Elastic Band) calculations using the pyiron framework with VASP. It includes
    classes and methods to manage initial and final states, reorder atoms,
    create NEB jobs, handle image interpolation, and collect and fit energy data.
    
    WARNING: It assumes remote and manual execution of VASP (on the cluster), it 
    only handles transfer of the input/output files, not the execution (as should 
    be ideally the case with pyiron).

Dependencies:
    - pyiron_atomistics
    - pyiron_base
    - NumPy
    - SciPy

Usage:
    from NEB import NEB
    # Create a NEB object, set initial and final states, and run the NEB calculation
    neb = NEB(pr, 'NEB_Calculation')
    neb.set_initial_state(initial_job)
    neb.set_final_state(final_job)
    neb_job = neb.create_neb_job(job_name='neb', nImages=5)
    neb_job.write_and_transfer_to_remote()
    # run the job manually on HPC
    neb_job.transfer_from_remote_and_collect()
    # get the barrier
    print(job['output/neb/barrier_forward'])

==============================================================================
"""

__author__ = "David Holec"
__version__ = "0.0.1"
__email__ = "david.holec@unileoben.ac.ay"
__status__ = "Prototype"  # Options: "Development", "Production", "Prototype"
__date__ = "2024-08-09"
__license__ = "BSD 3-Clause License"


from pyiron_atomistics.vasp.vasp import Vasp as Vasp_job
from pyiron_atomistics.vasp.structure import read_atoms
from pyiron_atomistics.vasp.parser.oszicar import Oszicar
from pyiron_atomistics.atomistics.job.atomistic import Trajectory
from pyiron_atomistics.atomistics.structure.atoms import Atoms
from pyiron_base import state
from pyiron_base.jobs.job.extension.jobstatus import job_status_successful_lst
from glob import glob
import os
import os.path
from shutil import rmtree
import numpy as np
from scipy.optimize import minimize
from numpy.polynomial import Polynomial

def get_species_symbols_new(self):
    """
    MONKEY PATCH to get correct order of POTCARs.
    
    Returns:
        numpy.ndarray: List of the symbols of the species.
    """
    sp = [self.elements[0]]
    for el in self.elements[1:]:
        if el != sp[-1]:
            sp.append(el)
    return np.array([el.Abbreviation for el in sp])

# Patch the `get_species_symbols` method of the `Atoms` class
Atoms.get_species_symbols = get_species_symbols_new


def _wrap_positions(struct, wrap):
    """
    Wrap atomic positions into the unit cell.

    Parameters:
        struct (Atoms): Atomic structure.
        wrap (float or list): Wrap value(s) for each dimension.

    Returns:
        Atoms: Structure with wrapped positions.
    """
    if isinstance(wrap, float):
        wrap = [wrap] * 3

    new_pos = []
    for pos in struct.get_scaled_positions():
        for i in range(3):
            if pos[i] > wrap[i]:
                pos[i] -= 1.0
        new_pos.append(pos)

    struct.set_scaled_positions(new_pos)
    return struct


class NEB():
    """
    Class to handle the setup and execution of NEB (Nudged Elastic Band) calculations.
    """
    
    def __init__(self, pr, NEB_calc_name):
        """
        Initialize the NEB class.

        Parameters:
            pr (Project): Pyiron project object.
            NEB_calc_name (str): Name of the NEB calculation.
        """
        self.pr = pr.create_group(NEB_calc_name)
        self.initial = None
        
    def add_initial_state(self, initial_job, delete_existing=False, wrap=0.98):
        """
        Add the initial state for the NEB calculation.

        Parameters:
            initial_job (Vasp_job): Initial state as a VASP job.
            delete_existing (bool, optional): Whether to delete existing files.
            wrap (float, optional): Wrap value for the atomic positions.
        
        Raises:
            TypeError: If `initial_job` is not of type `Vasp_job`.
        """
        if not isinstance(initial_job, Vasp_job):
            raise TypeError("Initial job must be a VASP job (pyiron_atomistics.vasp.vasp.Vasp)")
        
        self.initial_job = initial_job
        self.initial = initial_job.get_structure()
        if wrap is not False:
            self.initial = _wrap_positions(initial_job.get_structure(), wrap)
        
        # Setup symbolic links and directory structure
        rel_path_to_initial = os.path.relpath(self.initial_job.path, self.pr.path)
        cwd = os.getcwd()
        os.chdir(self.pr.path)
        if delete_existing:
            try:
                os.remove('initial.h5')
                rmtree('initial_hdf5')
            except:
                pass
        os.symlink(rel_path_to_initial+'.h5', 'initial.h5')
        os.mkdir('initial_hdf5')
        os.chdir('initial_hdf5')
        os.symlink(os.path.join('..', rel_path_to_initial+'_hdf5', initial_job.name), 'initial', target_is_directory=True)
        os.chdir(cwd)
        
    def set_initial_state(self, initial_job, wrap=0.98):
        """
        Set the initial state from an existing VASP job.

        Parameters:
            initial_job (Vasp_job): Initial state as a VASP job.
            wrap (float, optional): Wrap value for the atomic positions.
        """
        self.initial_job = initial_job
        struct = initial_job.get_structure()
        self.initial = struct
        if wrap is not False:
            self.initial = _wrap_positions(struct, wrap)
    
    def set_final_state(self, final_job, wrap=0.98):
        """
        Set the final state for the NEB calculation.

        Parameters:
            final_job (Vasp_job): Final state as a VASP job.
            wrap (float, optional): Wrap value for the atomic positions.
        """
        self.final_job = final_job
        struct = final_job.get_structure()
        self.final = struct
        if wrap is not False:
            self.final = _wrap_positions(struct, wrap)
        
        
    def reorder_final(self, moving_index, max_shift=0.05, verbose=True, wrap=0.98):
        """
        Reorder the final structure based on the initial structure.

        Parameters:
            moving_index (int): Index of the moving atom.
            max_shift (float, optional): Maximum allowable shift for atom mapping.
            verbose (bool, optional): If True, prints information about the atom positions.
            wrap (float, optional): Wrap value for the atomic positions.
        
        Raises:
            Exception: If there are issues with the atom mapping.
        """
        
        initial = self.initial
        final = self.final
        if wrap is not False:
            initial = _wrap_positions(initial, wrap)
            final = _wrap_positions(final, wrap)
            
        mapped = [False for _ in range(len(final))]
        mapping = [None for _ in range(len(initial))]
        
        for iat, at in enumerate(initial):
            if iat == moving_index and verbose:
                print('position in the initial structure:', at.position)
            else:
                nn = final.get_neighborhood(at.position, cutoff_radius=max_shift)
                if len(nn.distances) != 1:
                    raise Exception(f"Problem with mapping: atom {iat} in the initial structure was identified with {len(nn.distances)} images")
                elif mapped[nn.indices[0]]:
                    raise Exception(f"Problem with mapping: atom {iat} in the initial should be mapped to {nn.indices[0]} in the final structure but this atom has already been mapped.")
                else:
                    mapped[nn.indices[0]] = True
                    mapping[iat] = nn.indices[0]
        
        unmapped = [i for i in range(len(final)) if not mapped[i]]
        if len(unmapped) != 1:
            raise Exception(f"Problem with mapping: only 1 atom must remain after mapping! Unmapped atoms: {unmapped}")
        
        mapping[moving_index] = unmapped[0]
        if verbose:
            print('position in the final structure:', self.final.positions[unmapped[0]])
        
        new_positions = [final.positions[mapping[i]] for i in range(len(initial))]
        new_elements = [final.elements[mapping[i]] for i in range(len(initial))]
        new_final = self.pr.create.structure.atoms(
            cell=final.cell, 
            positions=new_positions, 
            elements=new_elements
        )
        self.final = new_final
        
    def create_neb_job(self, job_name='neb', nImages=3):
        """
        Create a NEB job.

        Parameters:
            job_name (str, optional): Name of the NEB job.
            nImages (int, optional): Number of images for the NEB calculation.
        
        Returns:
            NEBjob: Created NEB job object.
        """
        try:
            job = self.pr[job_name]
        except:
            job = self.pr.create_job(job_type=self.pr.job_type.Vasp, job_name=job_name)
            job.__class__ = NEBjob
            job.structure = self.initial         
            job.set_kpoints(mesh=[3, 3, 3], scheme='GC')
            job.calc_minimize(
                electronic_steps=60,
                retain_charge_density=False,
            )
            job.input.incar['ISIF'] = 2
            job.input.incar['ENCUT'] = 520  # from MP
            job.input.incar['IMAGES'] = nImages - 2
            job.input.incar['IBRION'] = 1
            job.input.incar['NFREE'] = 2
            job.save()
            
            job['input'].create_group('neb')            
            job['input/neb'].put('nImages', nImages)
            job['input/neb'].put('initial_image', self.initial)
            job['input/neb'].put('final_image', self.final)
            job['input/neb'].put('initial_energy', self.initial_job['output/generic/energy_pot'][-1])
            job['input/neb'].put('final_energy', self.final_job['output/generic/energy_pot'][-1])
            job._create_images()
            job.structure = job['input/neb/initial_image'].to_object()
            job.write_input()
        else:
            state.logger.warning(f'Job `{job_name}` exists. Loading from database instead of creating it!')
        return job


class NEBjob(Vasp_job):
    """
    Class for handling specific NEB (Nudged Elastic Band) job operations.
    """
    

    
    def _create_images(self, wrap=0.98):
        """
        Create interpolated images between the initial and final states.

        Parameters:
            wrap (float, optional): Wrap value for the atomic positions.
        """
        cells = []
        positions = []
        initial = _wrap_positions(self['input/neb/initial_image'].to_object(), wrap)
        final = _wrap_positions(self['input/neb/final_image'].to_object(), wrap)
        nImages = self['input/neb/nImages']
        
        for x in np.linspace(0, 1, nImages):
            cells.append(initial.cell * (1 - x) + final.cell * x)
            positions.append(initial.positions * (1 - x) + final.positions * x)
        
        self['input/neb'].put('initial_cells', cells)
        self['input/neb'].put('initial_positions', positions)
        
    def get_initial_images(self):
        """
        Retrieve the initial interpolated images.

        Returns:
            Trajectory: Trajectory object containing the initial images.
        """
        initial = self['input/neb/initial_image'].to_object()
        positions = self['input/neb/initial_positions']
        cells = self['input/neb/initial_cells']
        return Trajectory(positions, initial, cells=cells)
    
    def get_final_images(self):
        """
        Retrieve the final interpolated images after relaxation.

        Returns:
            Trajectory: Trajectory object containing the final images.
        """
        initial = self['input/neb/initial_image'].to_object()
        positions = self['output/neb/final_positions']
        cells = self['output/neb/final_cells']
        return Trajectory(positions, initial, cells=cells)
    
    def _write_images(self):
        """
        Write the interpolated images to the appropriate directories.
        """
        path = self.path + '_hdf5/' + self.job_name + '/'
        nImages = self['input/neb/nImages']
        initial_structs = self.get_initial_images()
        
        for i in range(nImages):
            os.chdir(path)
            case = f'{i:02d}'
            if not os.path.exists(case):
                os.mkdir(case)
            os.chdir(case)
            struct = initial_structs[i]
            struct.write('POSCAR', format='vasp')
    
    def write_and_transfer_to_remote(self):
        """
        Write the input files and transfer them to the remote cluster.

        Warnings:
            If the job is already marked as finished, a warning is issued and no action is taken.
        """
        if self.status not in job_status_successful_lst:
            self._write_images()
            self.write_input()
            
            filename = state.queue_adapter.convert_path_to_remote(
                path=self.project_hdf5.file_name
            )
            working_directory = state.queue_adapter.convert_path_to_remote(
                path=self.working_directory
            )
            for filename in glob(self.project_hdf5.path + '_hdf5/**', recursive=True):
                if os.path.isfile(filename):                
                    state.queue_adapter.transfer_file_to_remote(
                        file=filename, transfer_back=False
                    )
        else:
            state.logger.warning(f'Job `{self.job_name}` is marked as finished, uploading doesn\'t make sense')
            
    def collect_data(self, wrap=0.98):
        """
        Collect the data from the NEB calculation.

        Parameters:
            wrap (float or bool, optional): Wrap value for atomic positions. Set to `False` to disable wrapping.
        """
        path = self.path + '_hdf5/' + self.job_name + '/'
        nImages = self['input/neb/nImages']
        initial = self['input/neb/initial_image'].to_object()
        final = self['input/neb/final_image'].to_object()
        
        if wrap is not False:
            initial = _wrap_positions(initial, wrap)
            final = _wrap_positions(final, wrap)
        
        cells = [initial.cell]
        positions = [initial.positions]
        energy_pot = []
        
        for i in range(1, nImages - 1):
            case = f'{i:02d}'
            out = Oszicar()
            out.from_file(filename=path + case + '/OSZICAR')
            energy_pot.append(out.parse_dict['energy_pot'])
            struct = read_atoms(filename=path + case + '/CONTCAR')
            if wrap is not False:
                struct = _wrap_positions(struct, wrap)
            cells.append(struct.cell)
            positions.append(struct.positions)
        
        relaxation_steps = len(out.parse_dict['energy_pot'])
        energy_pot = np.array(
            [[self['input/neb/initial_energy']] * relaxation_steps] + 
            energy_pot +
            [[self['input/neb/final_energy']] * relaxation_steps]
        ).T
        
        cells.append(final.cell)
        positions.append(final.positions)
        
        self['output'].create_group('neb')
        self['output/neb'].put('energy_pot', energy_pot)
        self['output/neb'].put('final_cells', cells)
        self['output/neb'].put('final_positions', positions)
        ene = energy_pot[-1]
        self['output/neb'].put('barrier_forward', max(ene) - ene[0])
        self['output/neb'].put('barrier_backward', max(ene) - ene[-1])
        self._polynomial_fit_with_derivative_constraints(n=nImages + 1)
        
    def transfer_from_remote_and_collect(self):
        """
        Transfer files from the remote cluster and collect data locally.
        """
        self.transfer_from_remote()
        self.status.collect = True
        self.collect_data()
        self.status.finished = True
    
    def _polynomial_fit_with_derivative_constraints(self, n=4):
        """
        Fit a polynomial to the NEB energy profile with derivative constraints.

        The polynomial is constrained such that its derivatives at x=0 and x=1 are zero.

        Parameters:
            n (int, optional): The order of the polynomial.

        Returns:
            Polynomial: The fitted polynomial.
        """
        y = self['output/neb/energy_pot'][-1]
        x = np.linspace(0, 1, len(y))
        
        def poly_val_and_deriv(coeffs, x):
            """
            Evaluate the polynomial and its derivative.

            Parameters:
                coeffs (list or numpy.ndarray): Coefficients of the polynomial.
                x (float or numpy.ndarray): Point(s) at which to evaluate.

            Returns:
                tuple: Tuple containing the polynomial value and its derivative.
            """
            p = Polynomial(coeffs)
            dp = p.deriv()
            return p(x), dp(x)

        def objective(coeffs, x, y):
            """
            Objective function to minimize.

            Parameters:
                coeffs (list or numpy.ndarray): Coefficients of the polynomial.
                x (float or numpy.ndarray): x-values of the data points.
                y (float or numpy.ndarray): y-values of the data points.

            Returns:
                float: Sum of squared residuals.
            """
            p_vals, _ = poly_val_and_deriv(coeffs, x)
            return np.sum((y - p_vals) ** 2)

        def constraint_deriv_0(coeffs):
            """
            Constraint for the derivative at x=0 to be zero.

            Parameters:
                coeffs (list or numpy.ndarray): Coefficients of the polynomial.

            Returns:
                float: Derivative at x=0.
            """
            _, dp0 = poly_val_and_deriv(coeffs, 0)
            return dp0

        def constraint_deriv_1(coeffs):
            """
            Constraint for the derivative at x=1 to be zero.

            Parameters:
                coeffs (list or numpy.ndarray): Coefficients of the polynomial.

            Returns:
                float: Derivative at x=1.
            """
            _, dp1 = poly_val_and_deriv(coeffs, 1)
            return dp1

        def constraint_val_0(coeffs):
            """
            Constraint for the polynomial value at x=0 to match the first data point.

            Parameters:
                coeffs (list or numpy.ndarray): Coefficients of the polynomial.

            Returns:
                float: Difference between polynomial value at x=0 and y[0].
            """
            p0, _ = poly_val_and_deriv(coeffs, 0)
            return p0 - y[0]

        def constraint_val_1(coeffs):
            """
            Constraint for the polynomial value at x=1 to match the last data point.

            Parameters:
                coeffs (list or numpy.ndarray): Coefficients of the polynomial.

            Returns:
                float: Difference between polynomial value at x=1 and y[-1].
            """
            p1, _ = poly_val_and_deriv(coeffs, 1)
            return p1 - y[-1]

        initial_guess = np.ones(n + 1)

        def constraint_max_val(coeffs):
            """
            Constraint for the polynomial maximum to match the maximum y value.

            Parameters:
                coeffs (list or numpy.ndarray): Coefficients of the polynomial.

            Returns:
                float: Difference between maximum polynomial value and maximum y value.
            """
            p = Polynomial(coeffs)
            x_vals = np.linspace(0.2, 0.8, 10)
            p_vals = p(x_vals)
            return max(p_vals) - max(y)

        constraints = [
            {'type': 'eq', 'fun': constraint_deriv_0},
            {'type': 'eq', 'fun': constraint_deriv_1},
            {'type': 'eq', 'fun': constraint_val_0},
            {'type': 'eq', 'fun': constraint_val_1},
            {'type': 'eq', 'fun': constraint_max_val}
        ]

        result = minimize(objective, initial_guess, args=(x, y), constraints=constraints)
        optimized_coeffs = result.x
        fitted_polynomial = Polynomial(optimized_coeffs)
        self['output/neb'].put('fit', fitted_polynomial)
        
    def get_energy_fit(self):
        """
        Get the fitted polynomial for the NEB energy profile.

        Returns:
            Polynomial: The fitted polynomial.
        """
        return self['output/neb/fit']
