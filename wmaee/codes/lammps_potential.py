from typing import Optional, List, Tuple, Dict, Union
from wmaee.core.config import Config
import yaml
from os.path import join


def get_models() -> List[str]:
    """
    Get the list of interatomic potential models available in the LAMMPS potentials database.

    Returns:
        List[str]: A list of available potential models.
    """
    cfg = Config()
    models = cfg.get('applications').get('lammps').get('potentials')
    del models['root']
    return list(models.keys())



def get_potentials(model: str, species: Optional[List[str]] = None) -> Tuple[str, List[Dict[str, Union[str, List[str]]]]]:
    """
    Get LAMMPS potentials for a given model and optionally filtered by species.

    Parameters:
    - model (str): The name of the LAMMPS potential model.
    - species (Optional[List[str]]): A list of element symbols to filter the 
    potentials. If None, all potentials are returned.

    Returns:
    Tuple[str, List[Dict[str, Union[str, List[str]]]]]: A tuple containing the 
    potential root path and a list of potential dictionaries.
        - The first element of the tuple is the potential root path.
        - The second element is a list of potential dictionaries, where each 
        dictionary represents a LAMMPS potential and has the following keys:
            - 'elements' (List[str]): The list of element symbols associated 
            with the potential.
            - 'pot_file' (str): The filename of the potential file.
    """
    cfg = Config()
    root = cfg.get('applications').get('lammps').get('potentials').get('root')
    folder = cfg.get('applications').get('lammps').get('potentials').get(model)
    
    with open(join(root, folder, 'potentials.yaml')) as pot_file:
        potentials = dict(yaml.safe_load(pot_file))
    pots = []
    if species != None:
        for p in potentials['potentials']:
            if len(set(p['elements']) & set(species)) > 0:
                pots.append(p)
    else:
        # no species specified -> return all potentials
        pots = potentials['potentials']
    return potentials['pot_root'], pots