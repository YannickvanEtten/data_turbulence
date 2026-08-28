"""
1_download_hpc.py — stage-1 entrypoint: download ERA5 months for the CAT
pipeline on ADA, resumably and with a real integrity check.

This is the numbered stage entrypoint (like 2_diagnostics.py / 3_pipeline.py):
it DOES the download. The pure planning logic lives in download_plan.py, which
this imports. 1_download.py stays untouched as the original reference.

Guarantees (all from the Q-DOWNLOAD-1 audit, fail-loud throughout):
  - Never half-writes: each month downloads to <name>.tmp and is renamed to
    the final name only after it passes the integrity check. An interrupted
    run leaves a .tmp, never a truncated "final" file.
  - "Already done" means integrity-verified, not just os.path.exists. A month
    is skipped only if its final file exists AND re-passes the check.
  - Resumable/idempotent: a failed/interrupted month has no final file, so a
    re-run re-attempts it; a completed month is skipped. Safe to launch the
    same command (or SLURM array) repeatedly.
  - Fails loud: a downloaded file that fails the check is NOT renamed and
    raises; an existing final file that fails the check raises rather than
    being silently overwritten (use --force to deliberately re-download).

The integrity check verifies COUNTS (variables, pressure levels, timesteps)
against what the request actually asked for -- so the expected numbers can't
drift from the request. It does not yet check values; that's the diagnostics
layer's job.

NOTE (must be confirmed on the first real ADA grib): the dimension names used
to count levels/timesteps, and cfgrib's handling of all 7 variables in one
open, are asserted below but could not be exercised in the authoring
environment (no cfgrib/eccodes there). Run the trial-day mode interactively on
inter01-inter04 first and confirm the check passes on a known-good file before
trusting it at scale.
"""

import argparse
import os
from pathlib import Path

import xarray as xr

import download_plan


# ---------------------------------------------------------------------------
# CDS client -- imported lazily so the planning/verify functions (and their
# tests) don't require cdsapi to be installed.
# ---------------------------------------------------------------------------
def get_client(retry_max=8, sleep_max=60, timeout=120):
    """A CDS client that gives up in a bounded time.

    cdsapi's DEFAULTS ARE retry_max=500, sleep_max=120. On an unreachable CDS
    that is over sixteen hours of silent retrying -- longer than this job's
    entire walltime. A month that cannot be fetched would burn a compute node
    to the end of its allocation, print nothing useful, and be killed by the
    scheduler rather than reporting a cause.

    The bounded values here are the right shape for a SLURM array: 8 attempts
    backing off to a minute apart is roughly 8 minutes of genuine transient-
    error tolerance, then a loud failure. Because the downloader is resumable
    and idempotent -- a failed month leaves no final file, so re-submitting the
    same array retries exactly the months that failed -- failing fast is
    strictly better than hanging. Retrying is the job of the next submission,
    not of a wedged process.

    Discovered by accident, in a pre-flight probe that hung for five minutes
    against a blocked endpoint and would have kept going for sixteen hours.
    """
    import cdsapi
    return cdsapi.Client(retry_max=retry_max, sleep_max=sleep_max, timeout=timeout)


# ---------------------------------------------------------------------------
# Integrity check
# ---------------------------------------------------------------------------
# Candidate dimension names on a cfgrib-opened ERA5 pressure-level file. We
# search these rather than assume one, and RAISE if none is present, so a
# surprise layout fails loud instead of being silently miscounted.
LEVEL_DIM_CANDIDATES = ["isobaricInhPa", "pressure_level", "level", "plev"]
TIME_DIM_CANDIDATES = ["time", "valid_time", "step"]


def _find_dim(ds, candidates, what):
    for name in candidates:
        if name in ds.dims:
            return name
    raise KeyError(
        f"could not find the {what} dimension in {list(ds.dims)}; "
        f"looked for {candidates}. If ERA5's layout differs, add the real "
        f"name to the candidate list rather than guessing a default."
    )


def open_grib(path):
    """Open a GRIB file, failing loud on the common cfgrib multi-hypercube
    error rather than returning a silently-partial dataset."""
    try:
        # indexpath="" -> cfgrib keeps its index in memory and writes NO
        # <path>.<hash>.idx sidecar. We open each grib once, only to verify
        # counts, so a persisted index is pure clutter; and the download path
        # opens the .tmp then renames it, orphaning the sidecar. Suppressing at
        # the source means nothing is ever created OR deleted here, so this
        # cannot touch a real .grib.
        return xr.open_dataset(
            path, engine="cfgrib", backend_kwargs={"indexpath": ""}
        )
    except Exception as exc:
        raise RuntimeError(
            f"cfgrib could not open {path} as a single dataset ({exc}). "
            f"If this is the 'multiple values for key' error, the 7 variables "
            f"may need a filter_by_keys/backend_kwargs split -- resolve that "
            f"here explicitly; do not let it pass silently."
        ) from exc


def check_counts(ds, n_var, n_level, n_time):
    """Return a list of human-readable problems (empty list == passed).

    Counts, not names: robust to cfgrib short-names vs CDS long-names. Names
    are included in any problem string only to aid diagnosis.
    """
    problems = []

    got_var = len(ds.data_vars)
    if got_var != n_var:
        problems.append(
            f"variables: expected {n_var}, got {got_var} ({list(ds.data_vars)})"
        )

    level_dim = _find_dim(ds, LEVEL_DIM_CANDIDATES, "pressure-level")
    got_level = ds.sizes[level_dim]
    if got_level != n_level:
        problems.append(f"levels: expected {n_level}, got {got_level}")

    time_dim = _find_dim(ds, TIME_DIM_CANDIDATES, "time")
    got_time = ds.sizes[time_dim]
    if got_time != n_time:
        problems.append(f"timesteps: expected {n_time}, got {got_time}")

    return problems


def check_grid(ds, area):
    """Verify the horizontal grid matches the requested domain.

    THIS IS THE CHECK THAT DISTINGUISHES DOMAINS. A global month and a North
    Atlantic month have IDENTICAL variable, level and timestep counts -- 7, 3
    and 8*n_days for both -- so check_counts() passes on either. Only the grid
    size differs (36,421 points vs 1,038,240). Without this, a global file
    dropped into the North Atlantic series would verify clean.
    """
    problems = []
    exp_lat, exp_lon_options = download_plan.expected_grid(area)

    if "latitude" not in ds.sizes or "longitude" not in ds.sizes:
        problems.append(
            f"no latitude/longitude dimensions in {list(ds.sizes)}; cannot "
            f"confirm which domain this file covers"
        )
        return problems

    got_lat, got_lon = ds.sizes["latitude"], ds.sizes["longitude"]
    if got_lat != exp_lat:
        problems.append(f"latitudes: expected {exp_lat}, got {got_lat}")
    if got_lon not in exp_lon_options:
        expected = " or ".join(str(v) for v in sorted(exp_lon_options))
        problems.append(f"longitudes: expected {expected}, got {got_lon}")

    if problems:
        problems.append(
            f"-> grid mismatch means this file is probably for a DIFFERENT "
            f"DOMAIN than the one requested (area={area}). Check the filename's "
            f"domain code before overwriting anything."
        )
    return problems


def expected_counts(request):
    """Expected (variables, levels, timesteps) derived straight from the
    request that was downloaded -- so the check can't drift from the ask."""
    n_var = len(request["variable"])
    n_level = len(request["pressure_level"])
    n_time = download_plan.expected_timesteps(len(request["day"]))
    return n_var, n_level, n_time


def verify_file(path, request):
    """Open the file and return the list of integrity problems (empty == ok)."""
    ds = open_grib(path)
    try:
        n_var, n_level, n_time = expected_counts(request)
        problems = check_counts(ds, n_var, n_level, n_time)
        problems += check_grid(ds, request["area"])
        return problems
    finally:
        ds.close()


# ---------------------------------------------------------------------------
# Download one target (a month, or a trial day) to a final path, atomically.
# ---------------------------------------------------------------------------
def _fetch_verified(client, request, final, force=False):
    """Download `request` to `final`, atomically and integrity-checked.

    Returns "skip" if the final file already exists and re-passes the check,
    else "download". Raises loud on any integrity failure.
    """
    final = Path(final)
    tmp = final.with_name(final.name + ".tmp")

    if final.exists() and not force:
        problems = verify_file(final, request)
        if problems:
            raise RuntimeError(
                f"{final.name} already exists but FAILED integrity: {problems}. "
                f"Refusing to silently overwrite a supposedly-complete file. "
                f"Delete it by hand or re-run with --force to re-download."
            )
        return "skip"

    # Download to the temp name; the final name only appears after the check.
    client.retrieve(download_plan.DATASET, request, str(tmp))

    problems = verify_file(tmp, request)
    if problems:
        raise RuntimeError(
            f"downloaded {tmp.name} FAILED integrity: {problems}. Left the .tmp "
            f"in place for inspection; did NOT rename to {final.name}."
        )

    os.replace(tmp, final)  # atomic on the same filesystem
    return "download"


def download_month(client, year, month, out_dir, area=None, force=False,
                   domain=download_plan.DEFAULT_DOMAIN):
    request = download_plan.build_request(year, month, area=area, domain=domain)
    final = Path(out_dir) / download_plan.month_filename(year, month, domain)
    return _fetch_verified(client, request, final, force=force)


def subsampled_filename(year, month, days, domain):
    """Filename for a month downloaded on an explicit SUBSET of its days.

    Distinct from month_filename() for the same reason download_trial_day has
    its own name: a partial month must never be mistaken for a complete one by
    anything that globs the era5_<code>_YYYY-MM.grib series. The day list goes
    into the name rather than a sidecar, so the file is self-describing and a
    different sub-sample cannot silently overwrite this one.

        era5_glob_2000-01_d01-09-17-25.grib
    """
    code = download_plan.domain_code(domain)
    tag = "-".join(f"{int(d):02d}" for d in days)
    return f"era5_{code}_{year}-{month:02d}_d{tag}.grib"


def download_month_subset(client, year, month, days, out_dir, area=None,
                          force=False, domain=download_plan.DEFAULT_DOMAIN):
    """Download an explicit subset of a month's days.

    Used for the year-2000 global calibration pull, where thresholds are
    percentiles of a climatological distribution and therefore need coverage of
    the seasonal and diurnal cycles rather than every consecutive day. See
    jobs/02b_download_global_calib.sbatch.

    The integrity check needs no special-casing: expected_counts() derives the
    timestep count from len(request["day"]), so a 4-day request is checked
    against 4 x 8 = 32 timesteps, and check_grid() is unchanged.
    """
    request = download_plan.build_request(year, month, area=area, days=list(days),
                                          domain=domain)
    final = Path(out_dir) / subsampled_filename(year, month, days, domain)
    return _fetch_verified(client, request, final, force=force)


def download_trial_day(client, year, month, day, out_dir, area=None, force=False,
                       domain=download_plan.DEFAULT_DOMAIN):
    # Distinct filename so a one-day trial can never masquerade as a complete
    # month in the era5_<code>_YYYY-MM.grib series the pipeline treats as full
    # months. The domain code is included for the same reason it is on months.
    request = download_plan.build_request(year, month, area=area, days=[f"{day:02d}"],
                                          domain=domain)
    code = download_plan.domain_code(domain)
    final = Path(out_dir) / f"era5_trial_{code}_{year}-{month:02d}-{day:02d}.grib"
    return _fetch_verified(client, request, final, force=force)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Download ERA5 months for the CAT pipeline (resumable, "
                    "integrity-checked). Either a --start/--end span or a "
                    "single --trial-day."
    )
    p.add_argument("--out-dir", required=True,
                   help="Output directory, e.g. "
                        "/scistor/SBE-EDS-ClimateKoopman/yen230/era5")

    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--start", help="span start, 'YYYY-MM' (use with --end)")
    mode.add_argument("--trial-day", help="single trial day, 'YYYY-MM-DD'")
    p.add_argument("--end", help="span end, 'YYYY-MM' (required with --start)")

    p.add_argument("--index", type=int, default=None,
                   help="pick ONE month from the span by 0-based position "
                        "(for a SLURM array task); omit to process the whole span")
    p.add_argument("--domain", default=download_plan.DEFAULT_DOMAIN,
                   choices=sorted(download_plan.DOMAINS),
                   help="named domain to download. 'north_atlantic' is the trend "
                        "region; 'global' is the year-2000 calibration domain. The "
                        "domain code goes into every filename, so the two series "
                        "can share a directory without colliding.")
    p.add_argument("--area", type=float, nargs=4, default=None,
                   metavar=("N", "W", "S", "E"),
                   help="ad-hoc box override [N W S E]. Rarely needed -- prefer "
                        "--domain. If given it overrides the domain's box but NOT "
                        "the domain code in the filename, so use it only for "
                        "one-off experiments, never for production data.")
    p.add_argument("--days", default=None,
                   help="comma-separated day-of-month list, e.g. '1,9,17,25'. "
                        "Downloads only those days of each month in the span, "
                        "at full 3-hourly resolution, to a DIFFERENT filename "
                        "(era5_<code>_YYYY-MM_dNN-NN.grib) so a partial month "
                        "can never be mistaken for a complete one. Intended for "
                        "the year-2000 global calibration pull, where the "
                        "thresholds are percentiles of a climatological "
                        "distribution and need seasonal and diurnal coverage "
                        "rather than consecutive days. Days that a given month "
                        "does not have are dropped, so '1,9,17,25,31' is safe.")
    p.add_argument("--force", action="store_true",
                   help="re-download even if a final file already exists")

    args = p.parse_args(argv)
    if args.start and not args.end:
        p.error("--start requires --end")
    if args.days is not None:
        if args.trial_day:
            p.error("--days is for month spans; --trial-day is already one day")
        try:
            parsed = [int(tok) for tok in args.days.split(",") if tok.strip()]
        except ValueError:
            p.error(f"--days must be comma-separated integers, got {args.days!r}")
        if not parsed:
            p.error("--days was empty")
        bad = [d for d in parsed if not 1 <= d <= 31]
        if bad:
            p.error(f"--days values out of range 1..31: {bad}")
        args.days = sorted(set(parsed))
    return args


def main(argv=None):
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    area = args.area

    client = get_client()

    if args.trial_day:
        parts = args.trial_day.split("-")
        if len(parts) != 3:
            raise ValueError(f"expected --trial-day 'YYYY-MM-DD', got {args.trial_day!r}")
        year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
        result = download_trial_day(client, year, month, day, out_dir, area=area,
                                    force=args.force, domain=args.domain)
        print(f"trial {args.domain} {year}-{month:02d}-{day:02d}: {result}")
        return

    months = download_plan.month_span(args.start, args.end)
    if args.index is not None:
        if args.index < 0 or args.index >= len(months):
            raise IndexError(
                f"--index {args.index} out of range for span of {len(months)} "
                f"months ({args.start}..{args.end})"
            )
        months = [months[args.index]]

    for (year, month) in months:
        if args.days:
            # Drop days this month does not have, so one --days list works for
            # every month of the year (February included) without the job
            # script needing per-month logic.
            have = {int(d) for d in download_plan.month_days(year, month)}
            days = [d for d in args.days if d in have]
            if not days:
                raise ValueError(
                    f"none of --days {args.days} exist in {year}-{month:02d}"
                )
            if len(days) != len(args.days):
                dropped = sorted(set(args.days) - set(days))
                print(f"note: {year}-{month:02d} has no day(s) {dropped}; "
                      f"downloading {days}")
            result = download_month_subset(client, year, month, days, out_dir,
                                           area=area, force=args.force,
                                           domain=args.domain)
            print(f"{args.domain} {year}-{month:02d} days={days}: {result}")
        else:
            result = download_month(client, year, month, out_dir, area=area,
                                    force=args.force, domain=args.domain)
            print(f"{args.domain} {year}-{month:02d}: {result}")


if __name__ == "__main__":
    main()