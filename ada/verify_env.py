"""
ada/verify_env.py
=================
Environment checks, as a FILE rather than inline python.

WHY A FILE
----------
These checks were originally written as `pixi run python -c "...multi-line..."`
blocks inside the setup script. That fails: `pixi run` does not hand the string
straight to a shell, it parses it with its own task-shell parser first, and that
parser trips over ordinary English words that happen to be shell reserved words.
The first real run died on the word "in" -- inside a Python string, inside an
assertion message:

    Error: x failed to parse shell script ... Unsupported reserved word.
             in verification/ was produced at that rev.')

Nothing was wrong with the environment; the environment had just built cleanly.
A file has no quoting layer to get wrong, so this class of failure is gone.

Run:
    pixi run python ada/verify_env.py                 # all checks
    pixi run python ada/verify_env.py --skip-cds      # no network needed
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

EXPECTED_ROJAK_SUFFIX = "g25b8685c6"
REPO = Path(__file__).resolve().parent.parent

_failures: list[str] = []


def ok(msg: str) -> None:
    print(f"   [ok]   {msg}")


def bad(msg: str) -> None:
    print(f"   [FAIL] {msg}")
    _failures.append(msg)


def check_python() -> None:
    v = sys.version_info
    text = f"python {v.major}.{v.minor}.{v.micro}"
    if (v.major, v.minor) >= (3, 12):
        ok(text)
    else:
        bad(f"{text} -- rojak declares requires-python >=3.12")


def check_rojak(expected_suffix: str) -> None:
    try:
        import importlib.metadata as md
        version = md.version("rojak-cat")
    except Exception as exc:  # noqa: BLE001
        bad(f"rojak-cat not installed ({type(exc).__name__}: {exc})")
        return
    if version.endswith(expected_suffix):
        ok(f"rojak-cat {version} -- matches the pin in pixi.toml")
    else:
        bad(
            f"rojak-cat {version} does NOT end with {expected_suffix}. The "
            f"verification evidence under verification/ was produced at that "
            f"rev; a different one invalidates it."
        )


def check_grib_stack() -> None:
    try:
        import eccodes
        import cfgrib  # noqa: F401
        ok(f"cfgrib + ecCodes v{eccodes.codes_get_api_version()}")
    except Exception as exc:  # noqa: BLE001
        bad(f"GRIB stack unusable ({type(exc).__name__}: {exc})")


def check_diagnostics() -> None:
    """Load 2_diagnostics.py the way the pipeline does, and count the table."""
    try:
        spec = importlib.util.spec_from_file_location(
            "diagnostics", REPO / "2_diagnostics.py")
        module = importlib.util.module_from_spec(spec)
        # Registering BEFORE exec_module is mandatory: the module uses
        # `from __future__ import annotations`, so @dataclass resolves its
        # field annotations through sys.modules[cls.__module__].
        sys.modules["diagnostics"] = module
        spec.loader.exec_module(module)
    except Exception as exc:  # noqa: BLE001
        bad(f"2_diagnostics.py did not import ({type(exc).__name__}: {exc})")
        return
    n = len(module.REFERENCE_TABLE)
    if n == 21:
        ok(f"all {n} diagnostics import and register")
    else:
        bad(f"REFERENCE_TABLE has {n} entries, expected 21")


def check_cds() -> None:
    rc = Path.home() / ".cdsapirc"
    if not rc.exists():
        bad("~/.cdsapirc missing -- every download will fail")
        return
    mode = oct(rc.stat().st_mode)[-3:]
    url = ""
    for line in rc.read_text().splitlines():
        if line.lower().startswith("url"):
            url = line.strip()
            break
    if "api/v2" in url:
        bad(f"~/.cdsapirc uses the OLD CDS url ({url}). The new API needs "
            f"'url: https://cds.climate.copernicus.eu/api'")
        return
    ok(f"~/.cdsapirc present (mode {mode}), {url}")
    try:
        import cdsapi
        cdsapi.Client().check_authentication()
        ok("CDS authentication succeeded from this node")
    except Exception as exc:  # noqa: BLE001
        bad(f"CDS authentication failed ({type(exc).__name__}: {exc})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expect-rojak", default=EXPECTED_ROJAK_SUFFIX,
                        help="required suffix of the installed rojak version")
    parser.add_argument("--skip-cds", action="store_true",
                        help="skip the credential and authentication check")
    args = parser.parse_args()

    print(f"verifying the environment at {REPO}")
    check_python()
    check_rojak(args.expect_rojak)
    check_grib_stack()
    check_diagnostics()
    if args.skip_cds:
        print("   [skip] CDS checks (--skip-cds)")
    else:
        check_cds()

    print()
    if _failures:
        print(f"{len(_failures)} check(s) FAILED:")
        for f in _failures:
            print(f"  - {f}")
        return 1
    print("all environment checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
