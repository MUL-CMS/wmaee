from os import getcwd, mkdir, chdir
from os.path import exists

class working_directory(object):
    """
    Class for convenient change into a tree of subdirectories including their 
    creation.
    """
    
    def __split_path(self, path):
        """
        Helper function to split string into list of directories.
        """
        from os.path import split
        allparts = []
        while 1:
            parts = split(path)
            if parts[0] == path:  # sentinel for absolute paths
                allparts.insert(0, parts[0])
                break
            elif parts[1] == path: # sentinel for relative paths
                allparts.insert(0, parts[1])
                break
            else:
                path = parts[0]
                allparts.insert(0, parts[1])
        return allparts
    
    def __init__(self, path, create=True):
        """
        Parameters
        ----------
        path : str
            Name of the folder which we want to open for working.
        create : boolean, optional
            Specifies if the directory (tree) should be created in case it
            does not exist. The default is True.

        Returns
        -------
        None.

        """
        if not type(path) is str:
            raise TypeError('Path must be specified as a string')
        # split path into a list of directories
        self._path = self.__split_path(path)
        self._parent_dir = getcwd()
        self._create = create
        self._active = False

    def __enter__(self):        
        # we go through all subdirectories
        for sub in self._path:                        
            if not exists(sub):
                # if subdirectory doesn't exist, create it
                if self._create:
                    mkdir(sub)
                else:
                    raise Exception(f"The requested directory {sub} in {getcwd()}\
                                    doesn't exist and you do not want me to create it.")
            chdir(sub)
            self._active = True

    def __exit__(self, exc_type, exc_val, exc_tb):
        chdir(self._parent_dir)
        self._active = False

    @property
    def active(self):
        return self._active
    
# def grep(file, string):
#     with open(file, 'r') as fp:
#         for line in fp:
#             if string in line:
#                 yield line