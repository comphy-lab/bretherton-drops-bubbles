#!/usr/bin/env python3
"""
# Central-film observer and offline log acceptance

A Python port of `src-local/central-film.h`.  The in-solver observer is reset
on every restore, so a 72-hour continuation cannot inherit earlier windows.
This module rebuilds the same windows from one or more `c*-log` files.  Logs
that append across a compatible restart, or sequential logs from later
allocations, are stitched on physical time.  Overlapping or non-monotonic
prefix rows are skipped so a restart does not clear a hold that the film has
already earned.

This is a software contract for the observation predicate.  It does not
establish grid independence or agreement with Taylor's measurements.

## Author
Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
CoMPhy Lab, Durham University
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence, TextIO


@dataclass(frozen=True)
class CentralFilmConfig:
    advance_window: float
    minimum_observed_advance: float
    film_relative_tolerance: float
    front_speed_relative_tolerance: float
    minimum_coverage: float
    maximum_spatial_spread: float
    minimum_fresh_fraction: float
    minimum_cells: float
    shape_relative_tolerance: float
    hold_windows: int


@dataclass(frozen=True)
class CentralFilmMeasurement:
    time: float
    front: float
    rear: float
    centroid: float
    film: float
    spatial_spread: float
    coverage: float
    fresh_fraction: float
    minimum_cells: float


@dataclass
class CentralFilmWindow:
    advance: float = 0.0
    film: float = 0.0
    front_speed: float = 0.0
    rear_speed: float = 0.0
    centroid_speed: float = 0.0
    film_drift: float = 0.0
    front_speed_drift: float = 0.0
    shape_speed_mismatch: float = 0.0
    length_drift: float = 0.0
    minimum_coverage: float = 0.0
    maximum_spatial_spread: float = 0.0
    minimum_fresh_fraction: float = 0.0
    minimum_cells: float = 0.0
    hold_count: int = 0
    shape_hold_count: int = 0
    quality_ok: bool = False
    stable: bool = False
    film_converged: bool = False
    shape_steady: bool = False
    shape_converged: bool = False


@dataclass
class CentralFilmObserver:
    config: CentralFilmConfig
    first: CentralFilmMeasurement | None = None
    previous: CentralFilmMeasurement | None = None
    window_start: CentralFilmMeasurement | None = None
    film_integral: float = 0.0
    covered_advance: float = 0.0
    previous_window_film: float = 0.0
    previous_front_speed: float = 0.0
    minimum_coverage: float = 0.0
    maximum_spatial_spread: float = 0.0
    minimum_fresh_fraction: float = 0.0
    minimum_cells: float = 0.0
    hold_count: int = 0
    shape_hold_count: int = 0
    initialized: bool = False
    have_previous_window: bool = False


def central_film_quantile(values: Sequence[float], quantile: float) -> float:
    """Linearly interpolated quantile of a copied, sorted sample."""
    if (
        not values
        or not math.isfinite(quantile)
        or quantile < 0.0
        or quantile > 1.0
    ):
        return math.nan
    ordered = sorted(values)
    index = quantile * (len(ordered) - 1)
    lower = math.floor(index)
    lo = int(lower)
    hi = lo + 1 if lo + 1 < len(ordered) else lo
    return ordered[lo] + (index - lower) * (ordered[hi] - ordered[lo])


def central_film_median(values: Sequence[float]) -> float:
    return central_film_quantile(values, 0.5)


def central_film_observer_init(config: CentralFilmConfig) -> CentralFilmObserver:
    return CentralFilmObserver(config=config)


def central_film_observer_reset(observer: CentralFilmObserver) -> None:
    config = observer.config
    observer.__dict__.update(CentralFilmObserver(config=config).__dict__)


def central_film_relative_change(current: float, previous: float) -> float:
    return abs(current - previous) / abs(previous) if previous != 0.0 else math.inf


def _window_reset(observer: CentralFilmObserver, sample: CentralFilmMeasurement) -> None:
    observer.window_start = sample
    observer.previous = sample
    observer.film_integral = 0.0
    observer.covered_advance = 0.0
    observer.minimum_coverage = sample.coverage
    observer.maximum_spatial_spread = sample.spatial_spread
    observer.minimum_fresh_fraction = sample.fresh_fraction
    observer.minimum_cells = sample.minimum_cells


def _finite_sample(sample: CentralFilmMeasurement) -> bool:
    return all(
        math.isfinite(value)
        for value in (
            sample.time,
            sample.front,
            sample.rear,
            sample.centroid,
            sample.film,
            sample.spatial_spread,
            sample.coverage,
            sample.fresh_fraction,
            sample.minimum_cells,
        )
    )


def _sample_quality(config: CentralFilmConfig, sample: CentralFilmMeasurement) -> bool:
    return (
        sample.film > 0.0
        and sample.coverage >= config.minimum_coverage
        and sample.spatial_spread <= config.maximum_spatial_spread
        and sample.fresh_fraction >= config.minimum_fresh_fraction
        and (config.minimum_cells <= 0.0 or sample.minimum_cells >= config.minimum_cells)
    )


def central_film_observer_push(
    observer: CentralFilmObserver | None,
    sample: CentralFilmMeasurement,
    window: CentralFilmWindow | None,
) -> bool:
    """Add one sample. True only when an advance window has just closed."""
    if observer is None or window is None:
        return False
    if not _finite_sample(sample):
        central_film_observer_reset(observer)
        return False
    config = observer.config
    if not _sample_quality(config, sample):
        central_film_observer_reset(observer)
        return False
    if not observer.initialized:
        observer.initialized = True
        observer.first = sample
        _window_reset(observer, sample)
        return False
    assert observer.previous is not None and observer.window_start is not None
    if sample.time <= observer.previous.time or sample.front <= observer.previous.front:
        central_film_observer_reset(observer)
        observer.initialized = True
        observer.first = sample
        _window_reset(observer, sample)
        return False

    advance = sample.front - observer.previous.front
    if advance > 0.0:
        observer.film_integral += 0.5 * (sample.film + observer.previous.film) * advance
        observer.covered_advance += advance
    observer.minimum_coverage = min(observer.minimum_coverage, sample.coverage)
    observer.maximum_spatial_spread = max(
        observer.maximum_spatial_spread, sample.spatial_spread
    )
    observer.minimum_fresh_fraction = min(
        observer.minimum_fresh_fraction, sample.fresh_fraction
    )
    observer.minimum_cells = min(observer.minimum_cells, sample.minimum_cells)
    observer.previous = sample

    span = sample.front - observer.window_start.front
    elapsed = sample.time - observer.window_start.time
    if span < config.advance_window or observer.covered_advance <= 0.0 or elapsed <= 0.0:
        return False

    window.__dict__.update(CentralFilmWindow().__dict__)
    window.advance = span
    window.film = observer.film_integral / observer.covered_advance
    window.front_speed = span / elapsed
    window.rear_speed = (sample.rear - observer.window_start.rear) / elapsed
    window.centroid_speed = (
        sample.centroid - observer.window_start.centroid
    ) / elapsed
    window.minimum_coverage = observer.minimum_coverage
    window.maximum_spatial_spread = observer.maximum_spatial_spread
    window.minimum_fresh_fraction = observer.minimum_fresh_fraction
    window.minimum_cells = observer.minimum_cells

    length0 = observer.window_start.front - observer.window_start.rear
    length1 = sample.front - sample.rear
    window.length_drift = abs(length1 - length0) / length0 if length0 > 0.0 else math.inf
    speed_scale = abs(window.centroid_speed)
    window.shape_speed_mismatch = (
        max(
            abs(window.front_speed - window.centroid_speed),
            abs(window.rear_speed - window.centroid_speed),
        )
        / speed_scale
        if speed_scale > 0.0
        else math.inf
    )
    window.quality_ok = (
        window.film > 0.0
        and window.front_speed > 0.0
        and window.minimum_coverage >= config.minimum_coverage
        and window.maximum_spatial_spread <= config.maximum_spatial_spread
        and window.minimum_fresh_fraction >= config.minimum_fresh_fraction
        and (config.minimum_cells <= 0.0 or window.minimum_cells >= config.minimum_cells)
    )
    window.shape_steady = (
        window.shape_speed_mismatch <= config.shape_relative_tolerance
        and window.length_drift <= config.shape_relative_tolerance
    )

    if observer.have_previous_window:
        window.film_drift = central_film_relative_change(
            window.film, observer.previous_window_film
        )
        window.front_speed_drift = central_film_relative_change(
            window.front_speed, observer.previous_front_speed
        )
        assert observer.first is not None
        observed = sample.front - observer.first.front
        window.stable = (
            window.quality_ok
            and observed >= config.minimum_observed_advance
            and window.film_drift <= config.film_relative_tolerance
            and window.front_speed_drift <= config.front_speed_relative_tolerance
        )
        observer.hold_count = observer.hold_count + 1 if window.stable else 0
        observer.shape_hold_count = (
            observer.shape_hold_count + 1 if window.stable and window.shape_steady else 0
        )
    else:
        window.film_drift = math.inf
        window.front_speed_drift = math.inf
        observer.hold_count = 0
        observer.shape_hold_count = 0
        observer.have_previous_window = True

    observer.previous_window_film = window.film
    observer.previous_front_speed = window.front_speed
    window.hold_count = observer.hold_count
    window.shape_hold_count = observer.shape_hold_count
    window.film_converged = observer.hold_count >= config.hold_windows
    window.shape_converged = observer.shape_hold_count >= config.hold_windows
    _window_reset(observer, sample)
    return True


def parse_params(path: Path) -> dict[str, float | str]:
    """Read the runner-written `case.params` next to a log.

    The solver's only parameter path is `src-local/parse_params.h`. This
    function does not re-parse the source; it reads the copy the runner
    already wrote for that case.
    """
    values: dict[str, float | str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" not in line:
            continue
        key, text = line.split("=", 1)
        token = text.strip()
        try:
            values[key.strip()] = float(token)
        except ValueError:
            values[key.strip()] = token
    return values


def config_from_params(params: Mapping[str, float | str]) -> CentralFilmConfig:
    def number(name: str, default: float) -> float:
        value = params.get(name, default)
        return float(value)

    return CentralFilmConfig(
        advance_window=number("advWin", 0.25),
        minimum_observed_advance=number("advMin", 1.0),
        film_relative_tolerance=number("bTol", 2e-3),
        front_speed_relative_tolerance=number("speedTol", 0.02),
        minimum_coverage=number("filmCoverage", 0.95),
        maximum_spatial_spread=number("filmFlatTol", 0.05),
        minimum_fresh_fraction=number("freshFracMin", 0.98),
        minimum_cells=number("filmCells", 4.0),
        shape_relative_tolerance=number("shapeTol", 0.02),
        hold_windows=int(number("convHold", 3)),
    )


def _truthy(value: float | str, default: bool = True) -> bool:
    if isinstance(value, str):
        return value.strip().lower() not in {"0", "false", "no"}
    return bool(value) if value is not None else default


@dataclass(frozen=True)
class LogSample:
    path: Path
    row: int
    time: float
    front: float
    rear: float
    centroid: float
    film: float
    spatial_spread: float
    coverage: float
    fresh_fraction: float
    minimum_cells: float


def parse_log(path: Path) -> list[LogSample]:
    """Read one central-film case log. Older eight-column logs are rejected."""
    samples: list[LogSample] = []
    with path.open(encoding="utf-8") as handle:
        for row, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 15:
                raise ValueError(
                    f"{path}:{row} has {len(parts)} columns; "
                    "central-film logs need bCentral and minCells"
                )
            try:
                values = [float(item) for item in parts[:15]]
            except ValueError as error:
                raise ValueError(f"{path}:{row} is not numeric") from error
            samples.append(
                LogSample(
                    path=path,
                    row=row,
                    time=values[2],
                    front=values[5],
                    rear=values[6],
                    centroid=values[13],
                    film=values[8],
                    spatial_spread=values[9],
                    coverage=values[10],
                    fresh_fraction=values[11],
                    minimum_cells=values[12],
                )
            )
    if not samples:
        raise ValueError(f"no data rows in {path}")
    return samples


def _first_sample_time(path: Path) -> float:
    samples = parse_log(path)
    return samples[0].time


def stitch_logs(paths: Sequence[Path]) -> list[LogSample]:
    """Concatenate logs in order, dropping only each file's overlapping prefix.

    Later rows in the same file are passed through so a non-monotonic or
    non-finite sample can reset the observer, matching the C contract.
    """
    stitched: list[LogSample] = []
    for path in paths:
        in_prefix = True
        for sample in parse_log(path):
            if in_prefix and stitched and (
                sample.time <= stitched[-1].time
                or sample.front <= stitched[-1].front
            ):
                continue
            in_prefix = False
            stitched.append(sample)
    if not stitched:
        raise ValueError("stitched log is empty")
    return stitched


def _measurement(sample: LogSample) -> CentralFilmMeasurement:
    return CentralFilmMeasurement(
        time=sample.time,
        front=sample.front,
        rear=sample.rear,
        centroid=sample.centroid,
        film=sample.film,
        spatial_spread=sample.spatial_spread,
        coverage=sample.coverage,
        fresh_fraction=sample.fresh_fraction,
        minimum_cells=sample.minimum_cells,
    )


@dataclass
class Acceptance:
    film_converged: bool
    shape_converged: bool
    success: bool
    windows: int
    last_window: CentralFilmWindow | None
    last_time: float | None
    samples: int
    eligible_samples: int
    solver_tol: float


def accept_samples(
    samples: Sequence[LogSample],
    params: Mapping[str, float | str],
    *,
    require_shape: bool | None = None,
) -> Acceptance:
    """Replay stitched samples with a single initial burn, not a per-restart burn."""
    config = config_from_params(params)
    observer = central_film_observer_init(config)
    window = CentralFilmWindow()
    t_ramp = float(params.get("tRamp", 1.0))
    burn = float(params.get("restartBurnR", 0.0))
    require = (
        _truthy(params.get("requireShapeSteady", 1.0))
        if require_shape is None
        else require_shape
    )
    first_front: float | None = None
    windows = 0
    eligible = 0
    last: CentralFilmWindow | None = None
    last_time: float | None = None
    for sample in samples:
        if first_front is None:
            first_front = sample.front
        eligible_now = sample.time >= t_ramp and sample.front >= first_front + burn
        if not eligible_now:
            continue
        eligible += 1
        if central_film_observer_push(observer, _measurement(sample), window):
            windows += 1
            last = CentralFilmWindow(**window.__dict__)
            last_time = sample.time
    film_converged = observer.hold_count >= config.hold_windows
    shape_converged = observer.shape_hold_count >= config.hold_windows
    success = film_converged and (not require or shape_converged)
    return Acceptance(
        film_converged=film_converged,
        shape_converged=shape_converged,
        success=success,
        windows=windows,
        last_window=last,
        last_time=last_time,
        samples=len(samples),
        eligible_samples=eligible,
        solver_tol=float(params.get("solverTol", math.nan)),
    )


def discover_logs(case_dir: Path) -> list[Path]:
    logs = list(case_dir.glob("c*-log"))
    if not logs:
        raise FileNotFoundError(f"no c*-log in {case_dir}")
    return sorted(logs, key=lambda path: (_first_sample_time(path), path.as_posix()))


def accept_case(case_dir: Path, extra_logs: Sequence[Path] = ()) -> Acceptance:
    params_path = case_dir / "case.params"
    if not params_path.is_file():
        raise FileNotFoundError(f"no case.params in {case_dir}")
    logs = list(extra_logs) if extra_logs else discover_logs(case_dir)
    return accept_samples(stitch_logs(logs), parse_params(params_path))


def _window_json(window: CentralFilmWindow | None) -> dict[str, float | int | bool] | None:
    if window is None:
        return None
    payload = {}
    for key, value in window.__dict__.items():
        if isinstance(value, float) and not math.isfinite(value):
            payload[key] = None
        else:
            payload[key] = value
    return payload


def _print_acceptance(label: str, result: Acceptance, stream: TextIO) -> None:
    window = result.last_window
    film = f"{window.film:.6g}" if window else "none"
    speed = f"{window.front_speed:.6g}" if window else "none"
    stream.write(
        f"{label}: success={result.success} film={result.film_converged} "
        f"shape={result.shape_converged} windows={result.windows} "
        f"t={result.last_time} b={film} Ca_b={speed} "
        f"solverTol={result.solver_tol:.6g}\n"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Replay central-film windows from case logs without resetting on restart."
    )
    parser.add_argument("case_dirs", nargs="*", type=Path)
    parser.add_argument(
        "--logs",
        nargs="+",
        type=Path,
        help="explicit log files, stitched in this order",
    )
    parser.add_argument("--params", type=Path, help="case.params when using --logs")
    parser.add_argument("--json", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    reports: list[dict[str, object]] = []
    if args.logs:
        if args.params is None:
            parser.error("--logs requires --params")
        result = accept_samples(stitch_logs(args.logs), parse_params(args.params))
        reports.append({"label": "stitched", "result": result})
        if not args.json:
            _print_acceptance("stitched", result, sys.stdout)
    elif not args.case_dirs:
        parser.error("supply case directories or --logs")
    else:
        for case_dir in args.case_dirs:
            result = accept_case(case_dir)
            reports.append({"label": str(case_dir), "result": result})
            if not args.json:
                _print_acceptance(str(case_dir), result, sys.stdout)
    if args.json:
        payload = [
            {
                "label": item["label"],
                "success": item["result"].success,
                "film_converged": item["result"].film_converged,
                "shape_converged": item["result"].shape_converged,
                "windows": item["result"].windows,
                "last_time": item["result"].last_time,
                "samples": item["result"].samples,
                "eligible_samples": item["result"].eligible_samples,
                "solverTol": (
                    item["result"].solver_tol
                    if math.isfinite(item["result"].solver_tol)
                    else None
                ),
                "last_window": _window_json(item["result"].last_window),
            }
            for item in reports
        ]
        print(json.dumps(payload, indent=2))
    return 0 if all(item["result"].success for item in reports) else 1


if __name__ == "__main__":
    raise SystemExit(main())
