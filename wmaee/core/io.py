from os import getcwd, mkdir, chdir
from os.path import exists, split
from uuid import uuid4
from shutil import rmtree
from typing import Optional, Generator, List

class working_directory:
    """
    A convenience class which provides syntactic sugar, allowing the user to change directories.
    Can also be nested.
    
    Parameters
    ----------
    name : Optional[str], optional
        Name of the directory. If None is given, `os.getcwd()` will be used. Default is None.
    create : bool, optional
        Whether to create the directory if it doesn't exist. Default is True.
    delete : bool, optional
        Whether to delete the directory after a 'with' clause. Default is False.
    """
    def __init__(self, name: Optional[str] = None, create: bool = True, delete: bool = False):
        self._name: str = str(uuid4()) if not name else name
        self._delete: bool = delete
        self._create: bool = create
        self._curr_dir: str = getcwd()
        self._active: bool = False
        
    def _split_path(self, path: str) -> List[str]:
        """
        Helper function to split a string into a list of directories.
        
        Parameters
        ----------
        path : str
            The path to be split.

        Returns
        -------
        list[str]
            List of directories.
        """
        allparts: list[str] = []
        while 1:
            parts = split(path)
            if parts[0] == path:  # Sentinel for absolute paths
                allparts.insert(0, parts[0])
                break
            elif parts[1] == path:  # Sentinel for relative paths
                allparts.insert(0, parts[1])
                break
            else:
                path = parts[0]
                allparts.insert(0, parts[1])
        return allparts

    def __enter__(self) -> 'working_directory':
        """
        Enter the working directory.

        Returns
        -------
        working_directory
            The working_directory object.
        """
        path = self._split_path(self._name)
        for sub in path:
            if not exists(sub):
                if self._create:
                    mkdir(sub)
                else:
                    raise Exception(f"The requested directory {sub} in {getcwd()} "
                                    f"doesn't exist, and you do not want me to create it.")
            chdir(sub)
        self._active = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit the working directory.

        Parameters
        ----------
        exc_type : type
            Type of the exception.
        exc_val : any
            The exception value.
        exc_tb : traceback
            The traceback object.
        """
        chdir(self._curr_dir)
        if self._delete:
            rmtree(self._name)
        self._active = False

    @property
    def name(self) -> str:
        """
        Get the name of the working directory.

        Returns
        -------
        str
            The name of the working directory.
        """
        return self._name

    @property
    def active(self) -> bool:
        """
        Check if the working directory is active.

        Returns
        -------
        bool
            True if the working directory is active, False otherwise.
        """
        return self._active


# Utility functions
def grep(file: str, string: str) -> Generator[str, None, None]:
    """
    A simple Python implementation of the Linux `grep` command.
    It goes through a (text) file and returns all lines containing a desired string as a generator.

    Parameters
    ----------
    file : str
        Name of the file to be opened for reading and searching of the string.
    string : str
        String to be searched for.

    Yields
    ------
    str
        (Next) line containing the string `string` in the file `file`.
    """    
    with open(file, 'r') as f:
        for line in f:
            if string in line:
                yield line