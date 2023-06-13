
import os
import uuid
import shutil
import contextlib
from typing import Dict, Type

# this beautiful solution was taken from:
# https://stackoverflow.com/questions/2059482/temporarily-modify-the-current-processs-environment
@contextlib.contextmanager
def override_environ(*remove, **update):
    """
    Temporarily updates the ``os.environ`` dictionary in-place.

    The ``os.environ`` dictionary is updated in-place so that the modification
    is sure to work in all situations.

    :param remove: Environment variables to remove.
    :type remove: str
    :param update: Dictionary of environment variables and values to add/update.
    """
    env = os.environ
    update = update or {}
    remove = remove or []

    # List of environment variables being updated or removed.
    stomped = (set(update.keys()) | set(remove)) & set(env.keys())
    # Environment variables and values to restore on exit.
    update_after = {k: env[k] for k in stomped}
    # Environment variables and values to remove on exit.
    remove_after = frozenset(k for k in update if k not in env)

    try:
        env.update(update)
        [env.pop(k, None) for k in remove]
        yield
    finally:
        env.update(update_after)
        [env.pop(k) for k in remove_after]


class working_directory(object):
    """
    A convenience class which syntactic sugar, allowing the user to change the directories.
    Can also be nested.
    """

    def __init__(self, name=None, prefix=None, delete=False):
        """
        Constructs a working_directory object
        :param name: name of the directory if None is given `os.getcwd()` will be used (default is `None`)
        :type name: Optional[str]
        :param prefix: a prefix where to locate the directory (default is `None`)
        :type prefix: Optional[str]
        :param delete: whether to delete the directory after a with clause (default is `False`)
        :type delete: bool
        """
        self._name = str(uuid.uuid4()) if not name else name
        self._delete = delete
        self._curr_dir = os.getcwd()
        self._active = False
        if prefix is not None:
            self._name = os.path.join(prefix, self._name)

    def __enter__(self):
        if not os.path.exists(self._name):
            os.mkdir(self._name)
        os.chdir(self._name)
        self._active = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self._curr_dir)
        if self._delete:
            shutil.rmtree(self._name)
        self._active = False

    @property
    def name(self):
        return self._name

    @property
    def active(self):
        return self._active


def merge(*dicts: Dict[Any, Any], factory: Type = dict, **kwargs) -> Dict[Any, Any]:
    """
    Merge all specified dictionaries in {dicts} into a single one. Moreover, {kwargs} will be included. In case
    duplicate keys exist, the order of passing the dictionaries will determine which values will sustain for the keys.
    The last updated is carried out using {kwargs}

    :param dicts: the dictionaries to merge
    :type dicts: Dict[Any, Any]
    :param factory: the constructor of the mapping type to create (default is `dict`)
    :type factory: Type
    :return: a merged dictionary of type {factory}
    :rtype: Dict[Any, Any]
    """

    merged = factory()
    for dictionary in dicts:
        merged.update(dictionary)
    merged.update(kwargs)
    return merged


def _wrap_in_iterable(o, factory=tuple):
    return factory((o,))


def ensure_iterable(o, exclude=(str, bytes), factory=tuple):
    """
    wraps an object {o} into an iterable it is not and iterable. the type of the iterable is specified by {factory}
    :param o: the object to wrap
    :param exclude: type list of iterable objects which need wrapping (default is (str, bytes))
    :type exclude: Iterable[type]
    :param factory: the Iterable type in which {o} should be wrapped
    :type factory: type
    :rtype: Iterable
    """
    if isinstance(o, collections.abc.Iterable):
        return o if not isinstance(o, exclude) else _wrap_in_iterable(o, factory=factory)
    else:
        return _wrap_in_iterable(o, factory=factory)
