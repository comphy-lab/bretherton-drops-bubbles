"""Software-contract tests for runContinuation.py; they import no Basilisk."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import runContinuation as rc  # noqa: E402


class LadderTests(unittest.TestCase):
    def test_geometric_ladder_ascends_and_ends_on_stop(self):
        ladder = rc.geometric_ladder(0.05, 0.2, 1.5)
        self.assertAlmostEqual(ladder[0], 0.05)
        self.assertAlmostEqual(ladder[-1], 0.2)
        self.assertTrue(all(b > a for a, b in zip(ladder, ladder[1:])))

    def test_geometric_ladder_descends(self):
        ladder = rc.geometric_ladder(0.2, 0.05, 2.0)
        self.assertEqual(len(ladder), 3)
        self.assertAlmostEqual(ladder[1], 0.1)

    def test_fresh_station_has_no_frame_carry_over(self):
        base = {"MAXlevel": "8", "Uframe0": "0.3", "xTarget": "1"}
        p = rc.plan_station(base, case_no=1100, ca=0.05, seed=None,
                            maxlevel=None, renewals=8)
        self.assertEqual(p["CaPrev"], "0")
        self.assertNotIn("Uframe0", p)
        self.assertNotIn("xTarget", p)

    def test_continued_station_keeps_seed_frame_and_ramps_ca(self):
        seed = {"Ca_in": 0.05, "U": 0.068, "xTarget": 3.36, "t": 250.0,
                "length": 4.2}
        p = rc.plan_station({}, case_no=1101, ca=0.075, seed=seed,
                            maxlevel=9, renewals=8)
        self.assertEqual(float(p["CaPrev"]), 0.05)
        self.assertEqual(float(p["Uframe0"]), 0.068)
        self.assertEqual(p["MAXlevel"], "9")
        self.assertGreater(float(p["tmax"]), 250.0)

    def test_receipt_and_log_parsing(self):
        tmp = Path(self._testMethodName)
        tmp.mkdir(exist_ok=True)
        try:
            (tmp / "station.json").write_text(json.dumps({"status": "SUCCESS", "Ca_in": 0.05}))
            (tmp / "c1000-log").write_text("# CaseNo 1000\n# i dt t length dVol/Vol0 pExcess deltaTail minCells\n"
                                           "5 0.1 0.5 4.2 1e-4 3.7 0 4.1\n")
            self.assertEqual(rc.read_receipt(tmp)["status"], "SUCCESS")
            row = rc.last_log_row(tmp, 1000)
            self.assertEqual(float(row["length"]), 4.2)
        finally:
            for f in tmp.iterdir():
                f.unlink()
            tmp.rmdir()


if __name__ == "__main__":
    unittest.main()
