import os
import yaml
from typing import Optional, Any, Dict


def load_config(path: Optional[str] = None) -> Dict[str, Any]:
    """
    Loads a YAML config file from the specified path. If the path is not specified,
    the routine will try to load "~/.wmaee.conf.yaml" in case it exists. Otherwise,
    it uses the value from the environment variable `WMAEE_CONFIG_FILE`. If all options fail,
    a `FileNotFoundError` is raised.

    Parameters
    ----------
    path : Optional[str], optional
        Path to the config file.
    
    Returns
    -------
    Dict[str, Any]
        Configuration dictionary.
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
        return dict(yaml.safe_load(config_handle))


class Config:
    """
    Configuration class that loads settings from a YAML file.

    Parameters
    ----------
    path : Optional[str], optional
        Path to the config file.

    Attributes
    ----------
    _config : Dict[str, Any]
        Configuration dictionary.
    """
    def __init__(self, path: Optional[str] = None):
        self._config: Dict[str, Any] = load_config(path)

    def get(self, items: str) -> Any:
        """
        Retrieve a value from the configuration.

        Parameters
        ----------
        items : str
            Configuration item key.

        Returns
        -------
        Any
            Value associated with the specified key.
        """
        return self._config.get(items)