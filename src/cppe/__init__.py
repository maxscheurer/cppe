from . import cpp
from .cpp import *
from .cpp import __build_type__, __parallel__
from .pycppe.tensors import *
from .backend import get_backend, set_backend
from . import potfile as _potfile


def PotfileReader(*args, **kwargs):
    if get_backend() == "python":
        return _potfile.PotfileReader(*args, **kwargs)
    return cpp.PotfileReader(*args, **kwargs)


def NuclearFields(*args, **kwargs):
    if get_backend() == "python":
        from .fields import NuclearFields as _PythonNuclearFields

        return _PythonNuclearFields(*args, **kwargs)
    return cpp.NuclearFields(*args, **kwargs)


def BMatrix(*args, **kwargs):
    if get_backend() == "python":
        from .solver import BMatrix as _PythonBMatrix

        return _PythonBMatrix(*args, **kwargs)
    return cpp.BMatrix(*args, **kwargs)


def InducedMoments(*args, **kwargs):
    if get_backend() == "python":
        from .solver import InducedMoments as _PythonInducedMoments

        return _PythonInducedMoments(*args, **kwargs)
    return cpp.InducedMoments(*args, **kwargs)


def MultipoleExpansion(*args, **kwargs):
    if get_backend() == "python":
        from .multipole import MultipoleExpansion as _PythonMultipoleExpansion

        return _PythonMultipoleExpansion(*args, **kwargs)
    return cpp.MultipoleExpansion(*args, **kwargs)


def MultipoleFields(*args, **kwargs):
    if get_backend() == "python":
        from .multipole_fields import MultipoleFields as _PythonMultipoleFields

        return _PythonMultipoleFields(*args, **kwargs)
    return cpp.MultipoleFields(*args, **kwargs)


def CppeState(*args, **kwargs):
    if get_backend() == "python":
        from .state import CppeState as _PythonCppeState

        return _PythonCppeState(*args, **kwargs)
    return cpp.CppeState(*args, **kwargs)


all = [
    "__build_type__",
    "__parallel__",
    "cpp",
    "get_backend",
    "set_backend",
    "PotfileReader",
    "NuclearFields",
    "BMatrix",
    "InducedMoments",
    "MultipoleExpansion",
    "MultipoleFields",
    "CppeState",
]
