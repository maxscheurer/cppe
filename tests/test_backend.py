import importlib
from pathlib import Path

import pytest


def _reload_backend_module():
    import cppe.backend as backend

    return importlib.reload(backend)


def test_default_backend_is_cpp(monkeypatch):
    monkeypatch.delenv("CPPE_BACKEND", raising=False)
    backend = _reload_backend_module()
    assert backend.get_backend() == "cpp"


def test_backend_from_env(monkeypatch):
    monkeypatch.setenv("CPPE_BACKEND", "python")
    backend = _reload_backend_module()
    assert backend.get_backend() == "python"


def test_invalid_backend_env_falls_back_to_cpp(monkeypatch):
    monkeypatch.setenv("CPPE_BACKEND", "invalid")
    backend = _reload_backend_module()
    assert backend.get_backend() == "cpp"


def test_set_backend_changes_active_backend(monkeypatch):
    monkeypatch.delenv("CPPE_BACKEND", raising=False)
    backend = _reload_backend_module()

    backend.set_backend("python")
    assert backend.get_backend() == "python"

    backend.set_backend("cpp")
    assert backend.get_backend() == "cpp"


def test_set_backend_rejects_invalid_name(monkeypatch):
    monkeypatch.delenv("CPPE_BACKEND", raising=False)
    backend = _reload_backend_module()

    with pytest.raises(ValueError, match="Invalid backend"):
        backend.set_backend("fortran")


def test_cpp_namespace_is_available():
    import cppe

    assert hasattr(cppe, "cpp")
    assert hasattr(cppe.cpp, "PotfileReader")


def test_potfile_reader_dispatches_by_backend():
    import cppe
    import cppe.potfile as py_potfile

    potfile_path = Path(__file__).parent / "potfiles" / "pna_6w.pot"

    cppe.set_backend("cpp")
    cpp_reader = cppe.PotfileReader(str(potfile_path))
    assert isinstance(cpp_reader, cppe.cpp.PotfileReader)

    cppe.set_backend("python")
    py_reader = cppe.PotfileReader(str(potfile_path))
    assert isinstance(py_reader, py_potfile.PotfileReader)

    cppe.set_backend("cpp")
