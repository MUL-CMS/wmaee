import io
import os
import re
import contextlib
from frozendict import frozendict
from wmaee.core.utils import merge
from wmaee.units import HARTREE_TO_EV
from typing import Dict, Tuple, Any, TextIO, Union

ABINIT_PSEUDO_REGEX = re.compile(r"\"?(.*)[\",]?")
FLOAT_REGEX = r"[+-]?\d*\.\d+[eE]?[+-]?\d*"
INTEGER_REGEX = r"[+-]?\d+"
ABINIT_ENERGY_UNITS = '|'.join(("eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs"))


def identity(x):
    return x


def times(val: float):
    return lambda v: val*v


def parse_pseudo(path: str) -> Dict[str, str]:
    path, *_ = path.split(',')
    *_, xc, pps = os.path.basename(path).split('.')
    return dict(xc=xc, pps=pps)


def parse_ngkpt(*nums) -> Dict[str, Tuple[int, ...]]:
    return dict(ngkpt=tuple(map(int, nums)))


def parse_shiftk(*nums) -> Dict[str, Tuple[float, ...]]:
    return dict(shiftk=tuple(map(float, nums)))


from_hartree = times(HARTREE_TO_EV)
from_rydberg = times(HARTREE_TO_EV/2)
ENERGY_CONVERSION = dict(ev=identity, hartrees=from_hartree, hartree=from_hartree, ha=from_hartree, rydbergs=from_rydberg, rydberg=from_rydberg, ry=from_rydberg)

ABINIT_ENERGY_SCALAR_REGEX = re.compile(rf"({FLOAT_REGEX}|{INTEGER_REGEX})\s*({ABINIT_ENERGY_UNITS})?")


def make_parser(key: str, conversion, default: str):
    def parse_value(value: str, unit: str = default) -> Dict[str, float]:
        unit = unit.lower() if unit is not None else None

        try:
            value = float(value)
        except:
            raise ValueError(value)
        return {key: conversion.get(unit)(value) if unit is not None else value}

    return parse_value


ABINIT_KEYS = frozendict(
    pseudos=(parse_pseudo, ABINIT_PSEUDO_REGEX),
    ecut=(make_parser("ecut", ENERGY_CONVERSION, "eV"), ABINIT_ENERGY_SCALAR_REGEX),
    toldfe=(make_parser("toldfe", ENERGY_CONVERSION, "eV"), ABINIT_ENERGY_SCALAR_REGEX),
    nshiftk=(lambda v: dict(nshiftk=int(v.strip())), re.compile(f"({INTEGER_REGEX})")),
    kptopt=(lambda v: dict(kptopt=int(v.strip())), re.compile(f"({INTEGER_REGEX})")),
    nsppol=(lambda v: dict(nsppol=int(v.strip())), re.compile(f"({INTEGER_REGEX})")),
    ngkpt=(parse_ngkpt, re.compile(rf"({INTEGER_REGEX})\s+({INTEGER_REGEX})\s+({INTEGER_REGEX})")),
    shiftk=(parse_shiftk, re.compile(rf"({FLOAT_REGEX})\s+({FLOAT_REGEX})\s+({FLOAT_REGEX})"))
)


def parse_key(key: str, value: str):
    if key in ABINIT_KEYS:
        parser, regex = ABINIT_KEYS.get(key)
        return parser(*regex.match(value).groups())


def parse_file(handle):
    buf = io.StringIO()
    curr_key = None
    is_multiline_regex = re.compile(r"^([a-zA-Z]+)\s*$")
    is_command_regex = re.compile(r"^([a-zA-Z]+)")
    for line in handle:
        if line.startswith('#'):
            continue
        is_command = is_command_regex.match(line)
        is_multiline_command = is_multiline_regex.match(line)

        if is_command is None:
            buf.write(line)
        else:
            # check if there is something in the buffer
            if buf.tell() and buf.getvalue().strip():
                result = parse_key(curr_key, buf.getvalue())
                buf = io.StringIO()
                if result is not None:
                    yield result
            curr_key = is_command.group(0)
            if not is_multiline_command:
                _, value = line.split(' ', maxsplit=1)
                if '#' in value:  # strip off comments
                    value, _ = value.split('#', maxsplit=1)
                result = parse_key(curr_key, value)
                if result is not None:
                    yield result


def parse_abinit_input(path: Union[str, TextIO]) -> Dict[str, Any]:
    handle = open(path) if isinstance(path, str) else contextlib.nullcontext(path)
    with handle:
        return merge(*parse_file(handle))