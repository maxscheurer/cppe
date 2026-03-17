from .pycppe import *
from .pycppe import __build_type__, __parallel__
from .pycppe.tensors import *

try:
    from ._version import version as __version__
except Exception:
    try:
        from importlib.metadata import version as _pkg_version

        __version__ = _pkg_version("cppe")
    except Exception:
        __version__ = "0+unknown"

all = [
    "__build_type__",
    "__parallel__",
    "__version__",
]
