
import functools
from frozendict import frozendict


class UnmetRequirement(Exception):
    pass


class UnknownRequirement(Exception):
    pass


def test_gpaw() -> bool:
    try:
       from gpaw import GPAW
    except ImportError:
        return False
    else:
        return True


def test_kimpy() -> bool:
    try:
       import kimpy
    except ImportError:
        return False
    else:
        return True


def test_kim_query() -> bool:
    try:
       import kim_query
    except ImportError:
        return False
    else:
        return True


def test_nglview() -> bool:
    try:
       import nglview
    except ImportError:
        return False
    else:
        return True


REQUIREMENTS = frozendict(
    gpaw=test_gpaw,
    kimpy=test_kimpy,
    kim_query=test_kim_query,
    nglview=test_nglview
)


def requires(*requirements):
    """
    decorator that checks whether a "requirement" = module is avail able or not. Allows for optional dependencies
    in our package. Instead of an `ImportError` an `UnmetRequirement` exception is raised when the wrapped function
    is executed but the module is not available
    """

    def decorator(f):
        @functools.wraps(f)
        def test_requirements(*args, **kwargs):
            for req in requirements:
                if req not in REQUIREMENTS:
                    raise UnknownRequirement(req)
                elif not REQUIREMENTS[req]():
                    raise UnmetRequirement(req)
            return f(*args, **kwargs)

        return test_requirements
    return decorator