
from pymatgen.io.ase import AseAtomsAdaptor
from pyiron_base import PythonTemplateJob


class M3gnet(PythonTemplateJob):
    """
    M3GNet is a new materials graph neural network architecture that incorporates 3-body interactions.

    M3GNet (Materials Graph Neural Network) is designed to provide accurate interatomic potentials by 
    incorporating both 2-body and 3-body interactions, allowing for a more comprehensive representation 
    of the potential energy surface in materials science. It is particularly suited for simulations of 
    materials where complex interactions are significant.

    Literature:
    Chen, C., Ong, S.P. A universal graph deep learning interatomic potential for the periodic table. 
    Nat Comput Sci 2, 718–728 (2022). 
    https://doi.org/10.1038/s43588-022-00349-3
    """

    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        # Set default input parameters
        self.input.structure = None  # Input structure to be relaxed
        self.input.relax_cell = True  # Whether to relax the cell or only atomic positions
        self.input.fmax = 0.01  # Force tolerance for relaxation
        self.input.verbose = True  # Verbosity flag

        # we want to store the final structure for plotting reasons
        self.final_structure = None  # Initialize final_structure as an instance variable

    def format_numpy_array(self, numpy_list):
        """
        Converts a list of NumPy arrays to a multidimensional list.

        Parameters:
        - numpy_list (list of numpy.ndarray): List of NumPy arrays.

        Returns:
        - formatted_list (list of list of list): Multidimensional list.
        """
        numpy_array = np.array([f for f in numpy_list])
        formatted_list = numpy_array.tolist()
        return formatted_list

    def relax_structure_m3gnet(self):
        """
        Function to relax a given structure using a pre-trained MD potential.

        Returns:
        - final_structure (Structure): The relaxed structure.
        - relaxed_lattice_parameter (List): The relaxed lattice parameter.
        - energy_tot (List: flaot): The total energy for each step.

        """
        # Extract inputs
        structure = self.input.structure
        relax_cell = self.input.relax_cell
        fmax = self.input.fmax
        verbose = self.input.verbose

        # Initialize the relaxer with or without cell relaxation
        relaxer = Relaxer(relax_cell=relax_cell)

        # Perform relaxation
        relax_results = relaxer.relax(pyiron_to_ase(
            structure), fmax=fmax, verbose=verbose)

        # Extract the final relaxed structure
        final_structure = relax_results['final_structure']
        # convert to pyiron structure
        final_structure_pyiron = ase_to_pyiron(final_structure.to_ase_atoms())

        # Store the final structure as an instance variable
        self.final_structure = final_structure_pyiron

        # Calculate the relaxed lattice parameter
        relaxed_lattice_parameter = final_structure.lattice.abc

        # Get the full energy trajectory (list of energies)
        # Full list of energies
        energy_trajectory = relax_results['trajectory'].energies[:]

        # Get the full trajectory of the cell dimensions
        cells_trajectory = relax_results['trajectory'].cells
        formated_cells = self.format_numpy_array(cells_trajectory)

        # Get the full trajectory of the atomic positions
        positions_trajectory = relax_results['trajectory'].atom_positions
        formated_positions = self.format_numpy_array(positions_trajectory)

        # Get the full trajectory of the atomic positions
        force_trajectory = relax_results['trajectory'].forces
        formated_forces = self.format_numpy_array(force_trajectory)

        # Get the full trajectory of the stresses
        stress_trajectory = relax_results['trajectory'].stresses
        formated_stresses = self.format_numpy_array(stress_trajectory)

        # Print results
        if verbose:
            print(
                f"Relaxed lattice parameter is {relaxed_lattice_parameter} Å")
            print(f"Final energy is {energy_trajectory[-1]:.4f} eV")

        # Return results
        return relaxed_lattice_parameter, energy_trajectory, formated_positions, formated_cells, formated_forces, formated_stresses

    def run_static(self):
        """
        Runs the relaxation job and stores the results in the job output.
        """
        # Check if the structure is provided
        if self.input.structure is None:
            raise ValueError("Input structure not provided for relaxation.")

        # Run the relaxation using the relax_structure_m3gnet method
        relaxed_lattice_parameter, energy_trajectory, formated_positions, formated_cells, formated_forces, formated_stresses = self.relax_structure_m3gnet()

        # Store results in the output object
        self.output.cells = formated_cells
        self.output.positions = formated_positions
        self.output.relaxed_lattice_parameter = relaxed_lattice_parameter
        self.output.energy_tot = energy_trajectory
        self.output.forces = formated_forces
        self.output.stresses = formated_stresses

        # Mark the job as finished
        self.status.finished = True

        # Save the job to HDF5 file (persistent storage)
        self.to_hdf()

    def get_structure(self):
        """
        Retrieves the final relaxed structure after the job has been run.

        Returns:
        - final_structure (Structure): The relaxed structure.
        """
        if self.final_structure is None:
            raise ValueError(
                "The job has not been run or the structure is not available.")
        return self.final_structure


class Chgnet(PythonTemplateJob):
    """
    Crystal Hamiltonian Graph neural Network is pretrained on the GGA/GGA+U static and relaxation
    trajectories from Materials Project, a comprehensive dataset consisting of more than 1.5 Million
    structures from 146k compounds spanning the whole periodic table.

    CHGNet highlights its ability to study electron interactions and charge distribution in atomistic 
    modeling with near DFT accuracy. The charge inference is realized by regularizing the atom features 
    with DFT magnetic moments, which carry rich information about both local ionic environments and 
    charge distribution.

    Literature:
    Deng, Bowen and Zhong, Peichen and Jun, KyuJung and Riebesell, Janosh and Han, Kevin and Bartel, Christopher J. and Ceder, Gerbrand CHGNet as a pretrained universal neural network potential for charge-informed atomistic modelling 
    Nature Machine Intelligence, 1–11 (2023). 
    https://doi.org/10.1038/s42256-023-00716-3
    """

    def __init__(self, project, job_name):
        super().__init__(project, job_name)
        # Set default input parameters
        self.input.structure = None  # Input structure to be relaxed
        self.input.relax_cell = True  # Whether to relax the cell or only atomic positions
        self.input.fmax = 0.01  # Force tolerance for relaxation
        self.input.verbose = True  # Verbosity flag
        # "cpu", "cuda", or "mps", None selects automatically available option
        self.input.device = None
        self.input.steps = 500

        # we want to store the final structure for plotting reasons
        self.final_structure = None  # Initialize final_structure as an instance variable

    def format_numpy_array(self, numpy_list):
        """
        Converts a list of NumPy arrays to a multidimensional list.

        Parameters:
        - numpy_list (list of numpy.ndarray): List of NumPy arrays.

        Returns:
        - formatted_list (list of list of list): Multidimensional list.
        """
        numpy_array = np.array([f for f in numpy_list])
        formatted_list = numpy_array.tolist()
        return formatted_list

    def relax_structure_chgnet(self):
        """
        Function to relax a given structure using a pre-trained MD potential.

        Returns:
        - final_structure (Structure): The relaxed structure.
        - relaxed_lattice_parameter (List): The relaxed lattice parameter.
        - energy_tot (List: flaot): The total energy for each step.

        """
        # Extract inputs
        device = self.input.device
        structure = self.input.structure
        relax_cell = self.input.relax_cell
        fmax = self.input.fmax
        verbose = self.input.verbose
        steps = self.input.steps

        # Initialize the relaxer with or without cell relaxation
        relaxer = StructOptimizer(use_device=device)

        # Perform relaxation
        pymatgen_structure = AseAtomsAdaptor.get_structure(
            pyiron_to_ase(structure))
        relax_results = relaxer.relax(
            atoms=pymatgen_structure, verbose=verbose, fmax=fmax, steps=steps, relax_cell=relax_cell)

        # Extract the final relaxed structure
        final_structure = relax_results['final_structure']
        # convert to pyiron structure
        final_structure_pyiron = ase_to_pyiron(final_structure.to_ase_atoms())

        # Store the final structure as an instance variable
        self.final_structure = final_structure_pyiron

        # Calculate the relaxed lattice parameter
        relaxed_lattice_parameter = final_structure.lattice.abc

        # Get the full energy trajectory (list of energies)
        # Full list of energies
        energy_trajectory = relax_results['trajectory'].energies[:]

        # Get the full trajectory of the cell dimensions
        cells_trajectory = relax_results['trajectory'].cells
        formated_cells = self.format_numpy_array(cells_trajectory)

        # Get the full trajectory of the atomic positions
        positions_trajectory = relax_results['trajectory'].atom_positions
        formated_positions = self.format_numpy_array(positions_trajectory)

        # Get the full trajectory of the atomic positions
        force_trajectory = relax_results['trajectory'].forces
        formated_forces = self.format_numpy_array(force_trajectory)

        # Get the full trajectory of the stresses
        stress_trajectory = relax_results['trajectory'].stresses
        formated_stresses = self.format_numpy_array(stress_trajectory)

        # Get the full trajectory of the magmoms
        magmom_trajectory = relax_results['trajectory'].magmoms
        formated_magmoms = self.format_numpy_array(magmom_trajectory)

        # Print results
        if verbose:
            print(
                f"Relaxed lattice parameter is {relaxed_lattice_parameter} Å")
            print(f"Final energy is {energy_trajectory[-1]:.4f} eV")

        # Return results
        return relaxed_lattice_parameter, energy_trajectory, formated_positions, formated_cells, formated_forces, formated_stresses, formated_magmoms

    def run_static(self):
        """
        Runs the relaxation job and stores the results in the job output.
        """
        # Check if the structure is provided
        if self.input.structure is None:
            raise ValueError("Input structure not provided for relaxation.")

        # Run the relaxation using the relax_structure_m3gnet method
        relaxed_lattice_parameter, energy_trajectory, formated_positions, formated_cells, formated_forces, formated_stresses, formated_magmoms = self.relax_structure_chgnet()

        # Store results in the output object
        self.output.cells = formated_cells
        self.output.positions = formated_positions
        self.output.relaxed_lattice_parameter = relaxed_lattice_parameter
        self.output.energy_tot = energy_trajectory
        self.output.forces = formated_forces
        self.output.stresses = formated_stresses
        self.output.magmoms = formated_magmoms

        # Mark the job as finished
        self.status.finished = True

        # Save the job to HDF5 file (persistent storage)
        self.to_hdf()

    def get_structure(self):
        """
        Retrieves the final relaxed structure after the job has been run.

        Returns:
        - final_structure (Structure): The relaxed structure.
        """
        if self.final_structure is None:
            raise ValueError(
                "The job has not been run or the structure is not available.")
        return self.final_structure
