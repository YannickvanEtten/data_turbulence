"""Shared fixtures for the analytic verification suite."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))


def _load(module_name: str, filename: str):
    """Load a numbered module by path. Registering it in sys.modules BEFORE
    exec_module is mandatory: 2_diagnostics.py uses
    `from __future__ import annotations`, so @dataclass resolves its field
    annotations through sys.modules[cls.__module__]."""
    spec = importlib.util.spec_from_file_location(module_name, REPO / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def diag():
    return _load("diagnostics", "2_diagnostics.py")


@pytest.fixture(scope="session")
def grid():
    import synthetic as S
    return S.make_grid()


@pytest.fixture(scope="session")
def expect(grid):
    import synthetic as S
    return S.Expect(*grid)


@pytest.fixture(scope="session")
def prepared(diag):
    """The manufactured atmosphere, run through prepare_for_rojak."""
    import synthetic as S
    return diag.prepare_for_rojak(S.make_dataset())


@pytest.fixture(scope="session")
def rojak_out(diag, prepared):
    """All 14 rojak diagnostics, computed once for the whole session."""
    out, failures = diag.compute_rojak_diagnostics(prepared)
    assert not failures, f"rojak diagnostics failed: {[f.key for f in failures]}"
    return out
