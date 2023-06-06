import os
import logging
import shlex
import sys
from uuid import uuid4
from os import chdir, getcwd, mkdir
from os.path import join, exists
from shutil import rmtree
from wmaee.core.shell import Shell
from subprocess import PIPE, Popen
from wmaee.core.event import Event, EventHandler
from wmaee.utils import collection, unpack_single
from time import sleep, time as current_time
from io import TextIOWrapper, StringIO
from threading import Thread, Event as ThreadingEvent
from typing import Union, TextIO, Collection, Optional, Tuple, List, Dict, NoReturn, Iterator

def remove_white(string: str) -> str:
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


class REPL(Shell):

    def __init__(self, stdin=None, stdout=None, stderr=None, callbacks=None, forks=None, shell_cmd="/bin/bash"):
        """
        Create a shell object which holds the subprocess handle
        :param stdin: the input stream which is forwarded to handle.stdin
        :type stdin: InputStream
        :param stdout: the output stream where handle.stdout is forwarded to
        :type stdout: OutputStream
        :param stderr: the output stream where handle.stderr is forwarded to
        :type stderr: OutputStream
        :param callbacks: callbacks for the streams in order (in, out err)
        :type callbacks: Callbacks
        """
        super().__init__(stdin=stdin, stdout=stdout, stderr=stderr, callbacks=callbacks, forks=forks, shell_cmd=shell_cmd)
        self._magic = uuid4().hex


    async def execute_command(self, cmd, forks=None):
        self.handle.stdin.write(cmd.encode()+b'\n')
        await self.handle.stdin.drain()

    async def __call__(self, cmd):
        return await self.execute_command(cmd)