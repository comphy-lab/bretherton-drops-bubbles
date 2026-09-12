#!/usr/bin/env python3
"""Software tests for the Python central-film observer and log stitching."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "postProcess"))

from central_film import (  # noqa: E402
    CentralFilmConfig,
    CentralFilmMeasurement,
    CentralFilmWindow,
    accept_samples,
    central_film_median,
    central_film_observer_init,
    central_film_observer_push,
    central_film_quantile,
    discover_logs,
    parse_log,
    stitch_logs,
)


def standard_config() -> CentralFilmConfig:
    return CentralFilmConfig(
        advance_window=1.0,
        minimum_observed_advance=3.0,
        film_relative_tolerance=0.02,
        front_speed_relative_tolerance=0.02,
        minimum_coverage=0.8,
        maximum_spatial_spread=0.1,
        minimum_fresh_fraction=0.75,
        minimum_cells=4.0,
        shape_relative_tolerance=0.02,
        hold_windows=2,
    )


def measurement(time: float, front: float, film: float) -> CentralFilmMeasurement:
    return CentralFilmMeasurement(
        time=time,
        front=front,
        rear=front - 2.0,
        centroid=front - 1.0,
        film=film,
        spatial_spread=0.01,
        coverage=0.95,
        fresh_fraction=0.9,
        minimum_cells=6.0,
    )


def push(observer, window, time, front, film) -> bool:
    return central_film_observer_push(
        observer, measurement(time, front, film), window
    )


class QuantileTests(unittest.TestCase):
    def test_median_and_invalid_quantiles(self) -> None:
        self.assertEqual(central_film_median([9.0, 1.0, 5.0]), 5.0)
        self.assertEqual(central_film_median([4.0, 1.0, 3.0, 2.0]), 2.5)
        self.assertEqual(central_film_quantile([0.0, 10.0, 20.0, 30.0, 40.0], 0.25), 10.0)
        self.assertTrue(math.isnan(central_film_quantile([], 0.5)))
        self.assertTrue(math.isnan(central_film_quantile([1.0, 2.0], -0.01)))
        self.assertTrue(math.isnan(central_film_quantile([1.0, 2.0], 1.01)))
        self.assertTrue(math.isnan(central_film_quantile([1.0, 2.0], math.nan)))


class ObserverContractTests(unittest.TestCase):
    def test_stable_translating_film_converges(self) -> None:
        observer = central_film_observer_init(standard_config())
        window = CentralFilmWindow()
        self.assertFalse(push(observer, window, 0.0, 0.0, 0.1))
        self.assertTrue(push(observer, window, 1.0, 1.0, 0.1))
        self.assertTrue(window.quality_ok)
        self.assertFalse(window.stable)
        self.assertFalse(window.film_converged)
        self.assertTrue(push(observer, window, 2.0, 2.0, 0.1))
        self.assertEqual(window.hold_count, 0)
        self.assertTrue(push(observer, window, 3.0, 3.0, 0.1))
        self.assertTrue(window.stable)
        self.assertEqual(window.hold_count, 1)
        self.assertFalse(window.film_converged)
        self.assertTrue(push(observer, window, 4.0, 4.0, 0.1))
        self.assertEqual(window.hold_count, 2)
        self.assertTrue(window.film_converged)
        self.assertTrue(window.shape_steady)

    def test_relative_advance_after_a_late_anchor(self) -> None:
        observer = central_film_observer_init(standard_config())
        window = CentralFilmWindow()
        self.assertFalse(push(observer, window, 1000.0, 10000.0, 0.1))
        self.assertTrue(push(observer, window, 1001.0, 10001.0, 0.1))
        self.assertTrue(push(observer, window, 1002.0, 10002.0, 0.1))
        self.assertTrue(push(observer, window, 1003.0, 10003.0, 0.1))
        self.assertFalse(window.film_converged)
        self.assertTrue(push(observer, window, 1004.0, 10004.0, 0.1))
        self.assertTrue(window.film_converged)

    def test_quality_screens_clear_history(self) -> None:
        for field, value in (
            ("fresh_fraction", 0.2),
            ("coverage", 0.5),
            ("minimum_cells", 2.0),
            ("spatial_spread", 0.2),
            ("film", 0.0),
        ):
            observer = central_film_observer_init(standard_config())
            window = CentralFilmWindow()
            bad = CentralFilmMeasurement(
                **{**measurement(0.0, 0.0, 0.1).__dict__, field: value}
            )
            self.assertFalse(central_film_observer_push(observer, bad, window))
            self.assertFalse(observer.initialized)
            self.assertEqual(observer.hold_count, 0)

    def test_bad_quality_breaks_a_hold(self) -> None:
        observer = central_film_observer_init(standard_config())
        window = CentralFilmWindow()
        self.assertFalse(push(observer, window, 0.0, 0.0, 0.1))
        self.assertTrue(push(observer, window, 1.0, 1.0, 0.1))
        self.assertTrue(push(observer, window, 2.0, 2.0, 0.1))
        self.assertTrue(push(observer, window, 3.0, 3.0, 0.1))
        self.assertEqual(window.hold_count, 1)
        bad = CentralFilmMeasurement(**{**measurement(4.0, 4.0, 0.1).__dict__, "coverage": 0.5})
        self.assertFalse(central_film_observer_push(observer, bad, window))
        self.assertFalse(observer.initialized)
        self.assertFalse(push(observer, window, 5.0, 5.0, 0.1))
        self.assertTrue(push(observer, window, 6.0, 6.0, 0.1))
        self.assertFalse(window.film_converged)

    def test_drifting_film_is_not_stable(self) -> None:
        observer = central_film_observer_init(standard_config())
        window = CentralFilmWindow()
        self.assertFalse(push(observer, window, 0.0, 0.0, 0.1))
        self.assertTrue(push(observer, window, 1.0, 1.0, 0.1))
        self.assertTrue(push(observer, window, 2.0, 2.0, 0.2))
        self.assertGreater(window.film_drift, observer.config.film_relative_tolerance)
        self.assertFalse(window.stable)

    def test_film_and_shape_milestones_are_distinct(self) -> None:
        observer = central_film_observer_init(standard_config())
        window = CentralFilmWindow()
        for i in range(5):
            sample = measurement(float(i), float(i), 0.1)
            sample = CentralFilmMeasurement(
                **{
                    **sample.__dict__,
                    "rear": -2.0 + 0.5 * i,
                    "centroid": -1.0 + 0.8 * i,
                }
            )
            completed = central_film_observer_push(observer, sample, window)
            self.assertEqual(completed, i > 0)
        self.assertTrue(window.film_converged)
        self.assertFalse(window.shape_steady)
        self.assertFalse(window.shape_converged)

    def test_shape_convergence_needs_consecutive_windows(self) -> None:
        observer = central_film_observer_init(standard_config())
        window = CentralFilmWindow()
        for i in range(5):
            sample = CentralFilmMeasurement(
                **{
                    **measurement(float(i), float(i), 0.1).__dict__,
                    "rear": -2.0 + 0.5 * i,
                    "centroid": -1.0 + 0.8 * i,
                }
            )
            central_film_observer_push(observer, sample, window)
        self.assertTrue(window.film_converged)
        self.assertFalse(window.shape_steady)
        steady = CentralFilmMeasurement(
            **{**measurement(5.0, 5.0, 0.1).__dict__, "rear": 1.0, "centroid": 3.2}
        )
        self.assertTrue(central_film_observer_push(observer, steady, window))
        self.assertTrue(window.shape_steady)
        self.assertEqual(window.shape_hold_count, 1)
        self.assertFalse(window.shape_converged)
        steady = CentralFilmMeasurement(
            **{**measurement(6.0, 6.0, 0.1).__dict__, "rear": 2.0, "centroid": 4.2}
        )
        self.assertTrue(central_film_observer_push(observer, steady, window))
        self.assertTrue(window.shape_converged)

    def test_null_arguments_are_rejected(self) -> None:
        observer = central_film_observer_init(standard_config())
        window = CentralFilmWindow()
        sample = measurement(0.0, 0.0, 0.1)
        self.assertFalse(central_film_observer_push(None, sample, window))
        self.assertFalse(central_film_observer_push(observer, sample, None))


def _log_text(times: Sequence[float]) -> str:
    header = (
        "# CaseNo 1, MAXlevel 9, MINlevel 4, Ca 0.1, La 1, muR 0.01, "
        "rhoR 0.001, Rtube 0.7, Ldomain 16, solverTol 1e-6\n"
        "# i dt t ke dVol/Vol0 xTipF xTipR bFilm bCentral spatialSpread "
        "coverage freshFraction minCells centroid length\n"
    )
    rows = []
    for index, time in enumerate(times):
        front = time
        rear = time - 2.0
        centroid = time - 1.0
        rows.append(
            f"{index} 1e-4 {time:.6e} 1e-3 1e-6 {front:.6e} {rear:.6e} 0.1 "
            f"0.1 0.01 0.95 0.9 6.0 {centroid:.6e} 2.0\n"
        )
    return header + "".join(rows)


class LogStitchTests(unittest.TestCase):
    def test_overlapping_prefix_is_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            first = root / "c1-log"
            second = root / "c2-log"
            first.write_text(_log_text([0.0, 1.0, 2.0, 3.0]))
            second.write_text(_log_text([2.0, 3.0, 4.0, 5.0]))
            stitched = stitch_logs([first, second])
            self.assertEqual([sample.time for sample in stitched], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0])

    def test_short_second_allocation_does_not_converge_alone(self) -> None:
        params = {
            "advWin": 1.0,
            "advMin": 3.0,
            "bTol": 0.02,
            "speedTol": 0.02,
            "filmCoverage": 0.8,
            "filmFlatTol": 0.1,
            "freshFracMin": 0.75,
            "filmCells": 4.0,
            "shapeTol": 0.02,
            "convHold": 2,
            "tRamp": 0.0,
            "restartBurnR": 0.0,
            "requireShapeSteady": 1,
        }
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            first = root / "seg1-log"
            second = root / "seg2-log"
            first.write_text(_log_text([0.0, 1.0, 2.0, 3.0]))
            second.write_text(_log_text([3.0, 4.0]))
            first_only = accept_samples(parse_log(first), params)
            second_only = accept_samples(parse_log(second), params)
            stitched = accept_samples(stitch_logs([first, second]), params)
            self.assertFalse(first_only.success)
            self.assertFalse(second_only.success)
            self.assertFalse(second_only.film_converged)
            self.assertTrue(stitched.success)
            self.assertTrue(stitched.shape_converged)

    def test_logs_are_ordered_by_start_time_not_name(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "c10-log").write_text(_log_text([10.0, 11.0]))
            (root / "c2-log").write_text(_log_text([0.0, 1.0]))
            ordered = discover_logs(root)
            self.assertEqual([path.name for path in ordered], ["c2-log", "c10-log"])

    def test_only_the_overlapping_prefix_is_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            first = root / "c1-log"
            second = root / "c2-log"
            first.write_text(_log_text([0.0, 1.0, 2.0, 3.0]))
            second.write_text(_log_text([2.0, 3.0, 4.0, 3.5]))
            times = [sample.time for sample in stitch_logs([first, second])]
            self.assertEqual(times, [0.0, 1.0, 2.0, 3.0, 4.0, 3.5])

    def test_old_eight_column_log_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "c-old-log"
            path.write_text("# header\n0 1e-4 0.0 0.0 0.0 1.0 0.0 0.1\n")
            with self.assertRaises(ValueError):
                parse_log(path)


if __name__ == "__main__":
    unittest.main()
