
import os
import yaml
from jinja2 import Template
from typing import Optional, Any
from frozendict import frozendict


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


def singleton(class_, *args, **kwargs):
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

