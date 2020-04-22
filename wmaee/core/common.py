import os
import logging
from threading import Thread
from uuid import uuid4
from os import chdir, getcwd, mkdir
from os.path import join, exists
from shutil import rmtree


from typing import List

def remove_white(string):
    """
    Removes all whitespaces in a given string
    :param string: (str) the string
    :return: (str) a copy without whitespaces
    """
    whitespace = [' ', '\t', '\n']
    mystr = str(string)
    for removal in whitespace:
        mystr = mystr.replace(removal, '')
    return mystr


class ThreadWithReturnValue(Thread):
    """
    Small wraper around threading.Thread which stores the return value of the executed function
    """

    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args,
                                        **self._kwargs)

    def join(self, *args, **kwargs):
        Thread.join(self, *args, **kwargs)
        return self._return


class working_directory(object):
    """
    A convenience class which syntactic sugar, allowing the user to change the directories.
    Can also be nested.
    """

    def __init__(self, name=None, prefix=None, delete=False):
        """
        Constructs a working_directory object
        :param name: (str) name of the directory if None is given os.getcwd() will be used (default: None)
        :param prefix: (str) a prefix where to locate the directory (default: None)
        :param delete: (bool) wether to delete the directory after a with clause (default: False)
        """
        self._name = str(uuid4()) if not name else name
        self._delete = delete
        self._curr_dir = getcwd()
        self._active = False
        if prefix is not None:
            self._name = join(prefix, self._name)

    def __enter__(self):
        if not exists(self._name):
            mkdir(self._name)
        chdir(self._name)
        self._active = True

    def __exit__(self, exc_type, exc_val, exc_tb):
        chdir(self._curr_dir)
        if self._delete:
            rmtree(self._name)
        self._active = False

    @property
    def name(self):
        return self._name

    @property
    def active(self):
        return self._active


class LoggerMixin(object):
    """
    A mixin for logger
    """

    def __init__(self, **kwargs):
        super(LoggerMixin, self).__init__(**kwargs)

    @property
    def logger(self):
        return logging.getLogger(self.fullname())

    @classmethod
    def fullname(cls):
        """
        Returns the fully qualified name string of the cls
        :return: (str) the class identifier
        """
        name = '.'.join([
            cls.__module__,
            cls.__name__
        ])
        return name

def is_ipython():
    try:
        from IPython import get_ipython
    except ImportError:
        return False

    return get_ipython() is not None


if is_ipython():
    from tqdm import tqdm_notebook
    tqdm = tqdm_notebook
else:
    from tqdm import tqdm as tqdm_
    tqdm = tqdm_


def get_configuration_directory():
    """
    Build the path the to configuration directory for this module
    :return: (str) the absolute path of the configuration directory
    """
    if 'WMAEE_CONFIG_DIR' not in os.environ:
        if exists('.config'):
            return join(os.getcwd(), '.config')
        else:
            raise RuntimeError('No configuration directory found')
    else:
        return os.environ['WMAEE_CONFIG_DIR']