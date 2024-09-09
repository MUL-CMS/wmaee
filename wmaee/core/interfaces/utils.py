
import os
import yaml
import functools
from ase import Atoms
# from jinja2 import Template
from frozendict import frozendict
from wmaee.core.utils import ensure_iterable
from wmaee.core.requirements import requires
from typing import Optional, Any, List, Type, Callable, NoReturn, Union


def load_config(path: Optional[str] = None) -> frozendict:
    """
    Loads a YAML config file from path. If path is not specified, the routine will try to load "~/.wmaee.conf.yaml"
    in case it exists. Otherwise, the value from the environment variable `WMAEE_CONFIG_FILE`. If all options fail
    a `FileNotFoundError` is raised

    :param path: path to the config file
    :type path: Optional[str]
    :return: configuration dictionary
    :rtype: frozendict
    """
    if path is None:
        default_path = os.path.join(os.path.expanduser("~"), ".wmaee.conf.yaml")
        if os.path.exists(default_path):
            path = default_path
        else:
            env_path = os.environ.get("WMAEE_CONFIG_FILE")
            if env_path is None or not os.path.exists(env_path):
                raise FileNotFoundError("No config file was specified")
            else:
                path = env_path

    with open(path) as config_handle:
        return frozendict(yaml.safe_load(config_handle))


def singleton(class_: Type, *args, **kwargs) -> Callable[[], NoReturn]:
    """
    decorator to produce a singleton instance. `args` and `kwargs` are passed to the constructor to create the instance
    """
    instances = {}

    def getinstance():
        if class_ not in instances:
            instances[class_] = class_(*args, **kwargs)
        return instances[class_]

    return getinstance


class Config:

    def __init__(self, path: Optional[str] = None):
        self._config = load_config(path)

    def get(self, items: str) -> Any:
        return self._config.get(items)


# construct the configuration singleton object
Config = singleton(Config, path=None)


def render_command(application: str, **kwargs: Any) -> str:
    """
    Substitutes the parameters in `kwargs` to render a jinja2 template obtained from the configuration.

    :param application: name of the application to render the runscript for
    :type application: str
    :return: command string to launch the application
    :rtype: str
    """
    script = Config().get("applications").get(application).get("launch").get("script")
    args = Config().get("applications").get(application).get("launch").get("args")
    if not all(arg in kwargs or arg.lower() in kwargs for arg in args):
        raise ValueError(f"Missing at least one template argument: {args}")
    return Template(script).render({arg: kwargs.get(arg) or kwargs.get(arg.lower()) for arg in args})


@requires("kim_query")
def available_models(atoms: Atoms, species_logic: Union[str, List[str]] = "and", model_interface: Union[str, List[str]] = "any",
                     potential_type: Union[str, List[str]] = "any", simulator_name: Union[str, List[str]] = "any") -> List[str]:
    """
    Query the KIM online API for available potentials that might be used for {atoms}

    :param atoms: atoms objects
    :type atoms: ase.Atoms
    :param species_logic: reduction operation if {atoms} contains more than one species. Allowed choices are "and" and
        "or". (default is "and")
    :type species_logic: Union[str, List[str]]
    :param model_interface: model interface to select. Allowed choices are "sm" = Simulator Model, which refers to
        models bound to a specific MD code. "pm" = Portable Model(s) might be used with and MD code, and "any" if it
        does not matter (default is "any")
    :type model_interface: Union[str, List[str]]
    :param potential_type: the name or a list of names for potential types E.g. "meam", "eam" etc. (default is "any")
    :type potential_type: Union[str, List[str]]
    :param simulator_name: the name of a list of simulator names. E.g. "LAMMPS" (default is "any")
    :type simulator_name: Union[str, List[str]]
    :return: a list of available model names
    :rtype: List[str]
    """
    import kim_query
    species = list(set(atoms.symbols))
    model_interface = ensure_iterable(model_interface, factory=list)
    potential_type = ensure_iterable(potential_type, factory=list)
    simulator_name = ensure_iterable(simulator_name, factory=list)
    species_logic = ensure_iterable(species_logic, factory=list)
    return kim_query.get_available_models(species, species_logic=species_logic, model_interface=model_interface,
                                          simulator_name=simulator_name, potential_type=potential_type)


available_lammps_models = functools.partial(available_models, simulator_name="LAMMPS", model_interface="sm")
