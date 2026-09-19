#!/usr/bin/env python3
"""Ca continuation ladder for simulationCases/bretherton-comoving.c.

Each station is its own case directory (the case log pins one Ca per
directory). A station is seeded with a copy of the previous station's
``restart`` dump; the frame speed, centroid target and previous Ca travel in
``case.params`` (``Uframe0``, ``xTarget``, ``CaPrev``) because a dump holds
fields only. The inlet then ramps from ``CaPrev`` to ``Ca`` over ``tRamp`` and
the controller re-centres the bubble.

The ledger ``continuation.json`` carries one row per station with every field
of the station receipt plus route, seed and wall time. A run is resumed from
the ledger: a station directory holding a ``SUCCESS`` receipt is reused, never
re-run. A failed or incomplete station halves the step in log(Ca) from the last
converged station, at most ``--max-halvings`` times, then stops the ladder
with the reason recorded.

Uses only the standard library, so it runs on any controller.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "runSimulation.sh"


def read_params(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def write_params(path: Path, values: dict[str, str], header: str) -> None:
    lines = [f"# {header}"]
    lines += [f"{k}={v}" for k, v in values.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_ledger(path: Path) -> dict:
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"ladder": "Ca_in", "results": []}


def save_ledger(path: Path, ledger: dict) -> None:
    ledger["results"].sort(key=lambda r: (r.get("MAXlevel", 0), r["Ca_in"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    tmp.replace(path)


def read_receipt(case_dir: Path) -> dict | None:
    receipt = case_dir / "station.json"
    if not receipt.is_file():
        return None
    try:
        return json.loads(receipt.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def last_log_row(case_dir: Path, case_no: int) -> dict | None:
    log = case_dir / f"c{case_no}-log"
    if not log.is_file():
        return None
    header, last = None, None
    with log.open(encoding="utf-8") as fp:
        for line in fp:
            if line.startswith("# i "):
                header = line[2:].split()
            elif line and line[0].isdigit():
                last = line.split()
    if header is None or last is None or len(header) != len(last):
        return None
    return dict(zip(header, last))


def geometric_ladder(start: float, stop: float, factor: float) -> list[float]:
    if start <= 0 or stop <= 0 or factor <= 1:
        raise ValueError("start, stop must be positive and factor > 1")
    values = [start]
    if stop > start:
        while values[-1] * factor < stop * (1 - 1e-12):
            values.append(values[-1] * factor)
        values.append(stop)
    elif stop < start:
        while values[-1] / factor > stop * (1 + 1e-12):
            values.append(values[-1] / factor)
        values.append(stop)
    return values


def plan_station(base: dict[str, str], *, case_no: int, ca: float,
                 seed: dict | None, maxlevel: int | None,
                 renewals: float) -> dict[str, str]:
    """Parameters for one station given the seed receipt (or None)."""
    values = dict(base)
    values["CaseNo"] = str(case_no)
    values["Ca"] = f"{ca:.12g}"
    if maxlevel is not None:
        values["MAXlevel"] = str(maxlevel)
    if seed is None:
        values["CaPrev"] = "0"
        values.pop("Uframe0", None)
        values.pop("xTarget", None)
        return values
    # The dumped velocity field lives in the seed's frame, so the new station
    # must start with exactly the seed frame speed; the controller does the rest.
    values["CaPrev"] = f"{float(seed['Ca_in']):.12g}"
    values["Uframe0"] = f"{float(seed['U']):.17g}"
    values["xTarget"] = f"{float(seed['xTarget']):.17g}"
    # tmax is absolute on a restored clock: allow `renewals` film renewals of
    # the seed bubble at the predicted speed.
    length = float(seed.get("length", 0.0)) or 4.0
    u_pred = float(seed["U"]) * ca / float(seed["Ca_in"])
    values["tmax"] = f"{float(seed['t']) + renewals * length / u_pred:.6g}"
    return values


def run_station(case_dir: Path, params: Path, *, exec_name: str, threads: int,
                output_root: Path, dry_run: bool) -> int:
    cmd = ["bash", str(RUNNER), str(params), "--exec", exec_name,
           "--threads", str(threads)]
    env = dict(os.environ, OUTPUT_ROOT=str(output_root))
    print("RUN", " ".join(cmd), flush=True)
    if dry_run:
        return 0
    log = case_dir.parent / f"launch-{case_dir.name}.log"
    with log.open("a", encoding="utf-8") as fp:
        return subprocess.call(cmd, stdout=fp, stderr=subprocess.STDOUT, env=env)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", type=Path, required=True,
                        help="base parameter file for every station")
    parser.add_argument("--ca", type=float, nargs="*", default=None,
                        help="explicit Ca_in stations in ladder order")
    parser.add_argument("--ca-start", type=float)
    parser.add_argument("--ca-stop", type=float)
    parser.add_argument("--factor", type=float, default=1.5)
    parser.add_argument("--case-start", type=int, default=1100)
    parser.add_argument("--maxlevel", type=int, default=None)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--exec", dest="exec_name", default="bretherton-comoving.c")
    parser.add_argument("--output-root", type=Path, required=True,
                        help="registered run root; never the source tree")
    parser.add_argument("--ledger", type=Path, default=None)
    parser.add_argument("--seed-case", type=Path, default=None,
                        help="existing SUCCESS case directory to seed the first station")
    parser.add_argument("--renewals", type=float, default=8.0,
                        help="tmax budget per continued station in film renewals")
    parser.add_argument("--max-halvings", type=int, default=2)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.ca:
        stations = list(args.ca)
    elif args.ca_start and args.ca_stop:
        stations = geometric_ladder(args.ca_start, args.ca_stop, args.factor)
    else:
        parser.error("give --ca ... or --ca-start/--ca-stop")
    output_root = args.output_root.resolve()
    if ROOT in output_root.parents or output_root == ROOT:
        parser.error("output root must lie outside the source checkout")
    output_root.mkdir(parents=True, exist_ok=True)
    ledger_path = args.ledger or output_root / "continuation.json"
    ledger = load_ledger(ledger_path)
    base = read_params(args.base)
    for key in ("CaseNo", "Ca", "CaPrev", "Uframe0", "xTarget"):
        base.pop(key, None)

    seed: dict | None = None
    seed_dir: Path | None = None
    if args.seed_case:
        seed = read_receipt(args.seed_case)
        if not seed or seed.get("status") != "SUCCESS":
            parser.error(f"{args.seed_case} holds no SUCCESS receipt")
        row = last_log_row(args.seed_case, int(seed["CaseNo"]))
        if row:
            seed["length"] = float(row["length"])
        seed_dir = args.seed_case

    case_no = args.case_start
    pending = list(stations)
    halvings = 0
    exit_code = 0
    while pending:
        ca = pending[0]
        case_dir = output_root / str(case_no)
        receipt = read_receipt(case_dir)
        route = "fresh" if seed is None else "continuation"
        if receipt and receipt.get("status") == "SUCCESS" and \
                math.isclose(float(receipt["Ca_in"]), ca, rel_tol=1e-9):
            print(f"REUSE CaseNo {case_no} Ca_in={ca:g} (SUCCESS receipt)", flush=True)
        else:
            if case_dir.exists() and any(case_dir.iterdir()) and not receipt:
                print(f"STOP: {case_dir} is non-empty without a receipt; "
                      "inspect it before re-running", flush=True)
                return 2
            case_dir.mkdir(parents=True, exist_ok=True)
            if seed_dir is not None:
                shutil.copy2(seed_dir / "restart", case_dir / "restart")
            params = plan_station(base, case_no=case_no, ca=ca, seed=seed,
                                  maxlevel=args.maxlevel, renewals=args.renewals)
            params_path = output_root / f"station-{case_no}.params"
            write_params(params_path, params,
                         f"station CaseNo {case_no}: Ca_in={ca:g} route={route} "
                         f"seed={seed_dir.name if seed_dir else 'none'}")
            started = time.time()
            code = run_station(case_dir, params_path, exec_name=args.exec_name,
                               threads=args.threads, output_root=output_root,
                               dry_run=args.dry_run)
            wall = time.time() - started
            if args.dry_run:
                receipt = {"status": "DRY_RUN", "Ca_in": ca, "U": seed["U"] if seed else 0.0,
                           "xTarget": seed["xTarget"] if seed else 0.0, "t": 0.0,
                           "CaseNo": case_no}
            else:
                receipt = read_receipt(case_dir) or {"status": f"NO_RECEIPT_exit_{code}",
                                                     "Ca_in": ca, "CaseNo": case_no}
                receipt["wall_seconds"] = wall
        row = last_log_row(case_dir, case_no) if not args.dry_run else None
        if row:
            receipt["length"] = float(row["length"])
            receipt["dVol_over_Vol0"] = float(row["dVol/Vol0"])
            receipt["pExcess"] = float(row["pExcess"])
            receipt["deltaTail"] = float(row["deltaTail"])
            receipt["minCells"] = float(row["minCells"])
        receipt["route"] = route
        receipt["seed_case"] = seed_dir.name if seed_dir else None
        receipt["case_dir"] = str(case_dir)
        ledger["results"] = [r for r in ledger["results"]
                             if not (math.isclose(float(r["Ca_in"]), ca, rel_tol=1e-9)
                                     and r.get("MAXlevel") == receipt.get("MAXlevel"))]
        ledger["results"].append(receipt)
        save_ledger(ledger_path, ledger)
        print(f"STATION CaseNo {case_no} Ca_in={ca:g} status={receipt['status']} "
              f"Ca_b={receipt.get('Ca_b', float('nan')):.6g} "
              f"h/Rt={receipt.get('h_over_Rtube', float('nan')):.6f}", flush=True)

        ok = receipt["status"] in ("SUCCESS", "DRY_RUN")
        if ok:
            seed, seed_dir = receipt, case_dir
            pending.pop(0)
            halvings = 0
        else:
            if seed is None or halvings >= args.max_halvings:
                print(f"LADDER STOPS at Ca_in={ca:g}: {receipt['status']} "
                      f"after {halvings} halvings", flush=True)
                exit_code = 1
                break
            halvings += 1
            mid = math.sqrt(float(seed["Ca_in"]) * ca)
            print(f"HALVING step: inserting Ca_in={mid:g} before {ca:g}", flush=True)
            pending.insert(0, mid)
        case_no += 1
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
