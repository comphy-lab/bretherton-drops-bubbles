#!/usr/bin/env python3
"""Compare co-moving station receipts with the parent lab-frame film.

Reads one or more ``station.json`` receipts (or a ``continuation.json``
ledger) written by ``simulationCases/bretherton-comoving.c`` and prints the
deposited central film as ``h/R0`` and ``h/Rtube`` next to the parent
lab-frame values at the same finest cell size, plus the Aussillous-Quere
correlation evaluated at each code's own measured ``Ca_b``.

Standard library only. Figures are a separate task (publication-plots).
"""

from __future__ import annotations

import argparse
import csv
import re
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "postProcess" / "reference" / "parent-flat-film.csv"


def aussillous_quere(ca_b: float) -> float:
    """h/Rtube = 1.34 Ca^(2/3) / (1 + 2.5 * 1.34 Ca^(2/3))."""
    x = 1.34 * ca_b ** (2.0 / 3.0)
    return x / (1.0 + 2.5 * x)


def load_reference(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as fp:
        for row in csv.DictReader(line for line in fp if not line.startswith("#")):
            rows.append({k: (float(v) if k in ("level", "L0", "delta_over_R0", "Ca_in", "Ca_b", "H") else v)
                         for k, v in row.items()})
    return rows


def load_stations(paths: list[Path]) -> list[dict]:
    out = []
    for p in paths:
        # receipts written before the null-safe writer may contain bare nan
        data = json.loads(re.sub(r"\bnan\b", "null", p.read_text(encoding="utf-8")))
        if "results" in data:
            out.extend(r for r in data["results"] if r.get("status") == "SUCCESS")
        else:
            data.setdefault("case_dir", str(p.parent))
            out.append(data)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipts", nargs="+", type=Path,
                        help="station.json files or a continuation.json ledger")
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--json", type=Path, default=None,
                        help="write the comparison table as JSON")
    args = parser.parse_args(argv)

    reference = load_reference(args.reference)
    table = []
    for s in load_stations(args.receipts):
        delta = float(s["Ldomain"]) / 2 ** int(s["MAXlevel"])
        h = float(s["h"]); rt = float(s["Rtube"]); ca_b = float(s["Ca_b"])
        row = {
            "case": s.get("case_dir", ""), "status": s.get("status"),
            "Ca_in": float(s["Ca_in"]), "MAXlevel": int(s["MAXlevel"]),
            "delta_over_R0": delta, "Ca_b": ca_b,
            "h_over_R0": h, "h_over_Rtube": h / rt,
            "AQ_at_Ca_b": aussillous_quere(ca_b),
            "dev_from_AQ_percent": (100 * ((h / rt) / aussillous_quere(ca_b) - 1)
                                    if ca_b > 0 else float("nan")),
            "stagnant_film_residual": s.get("stagnant_film_residual"),
            "t_end": s.get("t"), "hold_count": s.get("hold_count"),
        }
        match = [r for r in reference if math.isclose(r["delta_over_R0"], delta, rel_tol=1e-6)
                 and math.isclose(r["Ca_in"], row["Ca_in"], rel_tol=1e-6)]
        if match:
            r = match[0]
            row.update({
                "parent_label": r["label"], "parent_Ca_b": r["Ca_b"],
                "parent_h_over_R0": r["H"] * 0.7, "parent_h_over_Rtube": r["H"],
                "parent_dev_from_AQ_percent": 100 * (r["H"] / aussillous_quere(r["Ca_b"]) - 1),
                "diff_h_percent": 100 * ((h / rt) / r["H"] - 1),
                "diff_Ca_b_percent": 100 * (ca_b / r["Ca_b"] - 1),
            })
        table.append(row)

    for row in table:
        print(f"Ca_in={row['Ca_in']:g} level={row['MAXlevel']} delta/R0={row['delta_over_R0']:.5g}")
        print(f"  co-moving : Ca_b={row['Ca_b']:.6f} h/R0={row['h_over_R0']:.6f} "
              f"h/Rt={row['h_over_Rtube']:.6f} AQ dev={row['dev_from_AQ_percent']:+.2f}% "
              f"stagnant={row['stagnant_film_residual']}")
        if "parent_label" in row:
            print(f"  parent {row['parent_label']}: Ca_b={row['parent_Ca_b']:.6f} "
                  f"h/R0={row['parent_h_over_R0']:.6f} h/Rt={row['parent_h_over_Rtube']:.6f} "
                  f"AQ dev={row['parent_dev_from_AQ_percent']:+.2f}%")
            print(f"  difference: h {row['diff_h_percent']:+.3f}%  Ca_b {row['diff_Ca_b_percent']:+.3f}%")
        else:
            print("  parent: no row at this delta and Ca_in in the reference table")
    if args.json:
        args.json.write_text(json.dumps(table, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
