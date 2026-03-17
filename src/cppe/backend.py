import os


_ALLOWED_BACKENDS = {"cpp", "python"}
_DEFAULT_BACKEND = "cpp"


def _normalize_backend(name: str) -> str:
    normalized = name.strip().lower()
    if normalized not in _ALLOWED_BACKENDS:
        allowed = ", ".join(sorted(_ALLOWED_BACKENDS))
        raise ValueError(f"Invalid backend '{name}'. Valid backends are: {allowed}.")
    return normalized


def _backend_from_env() -> str:
    raw = os.environ.get("CPPE_BACKEND", "")
    if not raw:
        return _DEFAULT_BACKEND
    try:
        return _normalize_backend(raw)
    except ValueError:
        return _DEFAULT_BACKEND


_backend = _backend_from_env()


def get_backend() -> str:
    return _backend


def set_backend(name: str) -> None:
    global _backend
    _backend = _normalize_backend(name)
