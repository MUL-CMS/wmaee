from pyiron_base import PythonTemplateJob
from pyiron_atomistics import Project, Atoms, pyiron_to_ase
from chgnet.model import CHGNet, StructOptimizer


class UPJob(PythonTemplateJob):
    def __init__(self, project, job_name):
        super(UPJob, self).__init__(project, job_name)
        self.input.atoms = None
        # its called fmax in chgnet and represents the force criteria for convergence
        self.input.ediffg = 0.01

    @property
    def structure(self):
        return self.input.atoms

    @structure.setter
    def structure(self, structure: Atoms):
        # Check if the structure passed is a pyiron `Atoms` object
        if isinstance(structure, Atoms):
            # if yes then convert to ase structure
            self.input.atoms = pyiron_to_ase(structure)
        elif isinstance(structure, ase.atoms.Atoms):
            # If the structure passed is already in ASE `Atoms` format, use it directly
            self.input.atoms = structure
        else:
            raise TypeError(
                f"Expected an pyiron or ase Atoms object but got {type(structure)}")

    @property
    def ediffg(self):
        return self.input.ediffg

    @ediffg.setter
    def ediffg(self, ediffg: float):
        self.input.ediffg = ediffg

    def run_static(self):
        print("running static calculation")
        relaxer = StructOptimizer()
        result = relaxer.relax(self.structure, fmax=self.ediffg, verbose=True)
        self.output.energy = self.structure.get_potential_energy()
        print("this is the output", self.structure.get_potential_energy())
        with self.project_hdf5.open("output/generic") as h5out:
            h5out["energy_pot"] = [self.output.energy]

        # Now we are finished
        self.status.finished = True
