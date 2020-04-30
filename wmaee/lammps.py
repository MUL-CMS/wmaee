
from wmaee.core.common import LoggerMixin

class LAMMPSInput(LoggerMixin):
    """
    A class used for creating lammps jobs
    """

    def view_potentials(cls, structure):
        """
        List all interatomic potentials for the current atomistic sturcture including all potential parameters.

        To quickly get only the names of the potentials you can use: self.potentials_list()

        Returns:
            pandas.Dataframe: Dataframe including all potential parameters.
        """
        from pyiron.lammps.potential import LammpsPotentialFile
        list_of_elements = set(structure.get_chemical_symbols())
        list_of_potentials = LammpsPotentialFile().find(list_of_elements)
        if list_of_potentials is not None:
            return list_of_potentials
        else:
            raise TypeError(
                "No potentials found for this kind of structure: ",
                str(list_of_elements),
            )


    def list_potentials(cls):
        """
        List of interatomic potentials suitable for the current atomic structure.

        use self.potentials_view() to get more details.

        Returns:
            list: potential names
        """
        return list(cls.view_potentials()["Name"].values)