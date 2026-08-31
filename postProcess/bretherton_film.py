#!/usr/bin/env python3
"""Extract film thickness and drop/bubble speed from bretherton.c case logs.

For each case directory the script reads ``c<CaseNo>-log`` (columns:
``i dt t ke dVol/Vol0 xTipF xTipR bFilm``), fits the front-tip velocity
over a trailing time window and averages the film thickness over the
same window. In code units (sigma = mu_c = R = 1) the tip velocity *is*
the capillary number ``Ca_b`` of the moving drop/bubble, so the measured
``b/Rtube`` can be compared directly with

- Bretherton (1961):        b/Rtube = 1.34 Ca_b^(2/3)
- Aussillous & Quere (2000): b/Rtube = 1.34 Ca_b^(2/3) / (1 + 3.35 Ca_b^(2/3))

The script writes one CSV row per case and prints a human-readable
summary. It performs no plotting.

Usage:
    python3 postProcess/bretherton_film.py simulationCases/1000 [more dirs...]
    python3 postProcess/bretherton_film.py --window 0.25 --out film.csv dirs...
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def read_log(case_dir: Path):
    """Return (header_params, rows) from the case log."""
    logs = sorted(case_dir.glob("c*-log"))
    if not logs:
        raise FileNotFoundError(f"no c<CaseNo>-log in {case_dir}")
    log = logs[0]
    params = {}
    rows = []
    with log.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                for key, value in re.findall(r"(\w+)\s+([-+0-9.eE]+)", line):
                    params[key] = float(value)
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    rows.append([float(v) for v in parts[:8]])
                except ValueError:
                    continue
    if not rows:
        raise ValueError(f"no data rows in {log}")
    return params, rows


def linear_fit(xs, ys):
    """Least-squares slope and intercept."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0.0:
        return 0.0, my
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return slope, my - slope * mx


def analyse_case(case_dir: Path, window: float):
    params, rows = read_log(case_dir)
    if "Rtube" not in params:
        raise ValueError("no Rtube in the log header; cannot normalise the film")
    rtube = params["Rtube"]
    ca_imposed = params.get("Ca", float("nan"))

    t_all = [r[2] for r in rows]
    t_end = t_all[-1]
    t_start = t_end * (1.0 - window)
    tail = [r for r in rows if r[2] >= t_start]
    if len(tail) < 5:
        tail = rows[-max(5, len(rows) // 4):]

    ts = [r[2] for r in tail]
    xtips = [r[5] for r in tail]
    bfilms = [r[7] for r in tail]

    ca_b, _ = linear_fit(ts, xtips)  # tip velocity = Ca_b in code units
    b_mean = sum(bfilms) / len(bfilms)
    b_over_r = b_mean / rtube

    bret = 1.34 * ca_b ** (2.0 / 3.0) if ca_b > 0 else float("nan")
    aq = bret / (1.0 + 2.5 * bret) if ca_b > 0 else float("nan")

    return {
        "case": case_dir.name,
        "Ca_imposed": ca_imposed,
        "Ca_b": ca_b,
        "Rtube": rtube,
        "b": b_mean,
        "b_over_Rtube": b_over_r,
        "bretherton": bret,
        "aussillous_quere": aq,
        "dev_bretherton": (b_over_r - bret) / bret if bret and bret > 0 else float("nan"),
        "dev_aq": (b_over_r - aq) / aq if aq and aq > 0 else float("nan"),
        "t_window": f"{t_start:.3g}..{t_end:.3g}",
        "n_samples": len(tail),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case_dirs", nargs="+", type=Path,
                        help="case directories, e.g. simulationCases/1000")
    parser.add_argument("--window", type=float, default=0.25,
                        help="trailing time-window fraction for the fit "
                             "(default: 0.25)")
    parser.add_argument("--out", type=Path, default=None,
                        help="optional CSV output path")
    args = parser.parse_args(argv)

    if not 0.0 < args.window < 1.0:
        parser.error("--window must be in (0, 1)")

    results = []
    for case_dir in args.case_dirs:
        try:
            results.append(analyse_case(case_dir, args.window))
        except (FileNotFoundError, ValueError) as exc:
            print(f"SKIP {case_dir}: {exc}", file=sys.stderr)

    if not results:
        print("No cases analysed.", file=sys.stderr)
        return 1

    fields = list(results[0].keys())
    if args.out:
        try:
            with args.out.open("w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=fields)
                writer.writeheader()
                writer.writerows(results)
        except OSError as exc:
            print(f"ERROR: cannot write {args.out}: {exc}", file=sys.stderr)
            return 1
        print(f"Wrote {args.out}")

    for r in results:
        print(f"{r['case']}: Ca_b = {r['Ca_b']:.4g}, "
              f"b/Rtube = {r['b_over_Rtube']:.4g} "
              f"(Bretherton {r['bretherton']:.4g}, "
              f"AQ {r['aussillous_quere']:.4g}; "
              f"dev {100 * r['dev_bretherton']:+.1f}% / "
              f"{100 * r['dev_aq']:+.1f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
