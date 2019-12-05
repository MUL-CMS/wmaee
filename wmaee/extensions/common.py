
from threading import Thread
from io import StringIO
from uuid import uuid4
from time import sleep
from os import chdir, getcwd
from shutil import rmtree
from os.path import join, exists
from os import mkdir
import logging

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

class StringStream(StringIO):
    """
    A class representing a dummy stream, which can be used to write data to and read from it
    """

    def __init__(self, string=''):
        """
        Create a string stream with a initial value
        :param string: (str) the initial value (default: "")
        """
        super(StringStream, self).__init__(initial_value=string)
        self._pos = 0
        self._remaining = 0
        self._length = 0

    def read(self, size=-1):
        """
        Performs a read operation on the StringStream object. Blocks if not enough data is available
        :param size: (int) number of characters to be read from the stream (default: -1)
        :return: (str) data read from the StringStream object
        """
        while self._remaining < size:
            sleep(0.01)
        result = super(StringStream, self).read()
        # Increase position, from current position seek( ..., 1)
        result_length = len(result)
        # Increase position, and cosume
        self._pos += result_length
        self._remaining = self._length - self._pos
        self.seek(self._pos)
        return result

    def write(self, s):
        """
        Write data to the StringStream object
        :param s: (str) the data
        """
        write_length = len(s)
        super(StringStream, self).write(s)
        # After write file is at the end
        # Seek from back and make it available
        self._length += write_length
        self._remaining = self._length - self._pos

    def readline(self, size=-1, block=True):
        """
        Reads a line from the StringStream object. Block if not a full line is available
        :param size: (int) number of characters to read (default: -1)
        :param block: (bool) wether to block until  a line is available (default: True)
        :return: (str) the data read from the StringStream object
        """
        if self.tell() != self._pos:
            self.seek(self._pos)
        result = super(StringStream, self).readline()
        result_length = len(result)

        self._pos += result_length
        self._remaining = self._length - self._pos
        # Seek new position
        self.seek(self._pos)
        # if block:
        #    while not result:
        #        result = super(StringStream, self).readline()
        #        sleep(0.025)
        return result


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