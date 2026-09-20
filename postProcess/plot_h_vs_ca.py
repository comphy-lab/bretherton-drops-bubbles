#!/usr/bin/env python3
"""Deposited film h/R0 against bubble Ca_b for the co-moving ladders.

Reads one or more ``continuation.json`` ledgers written by ``runContinuation.py``
(one per cell size), keeps SUCCESS rows, and plots them log-log with the
Aussillous-Quere and Bretherton laws, the parent lab-frame points and the FEM
sweep. The rendered figure carries only axis labels and legend entries;
interpretation belongs in the caption of the embedding document.

    python3 postProcess/plot_h_vs_ca.py --ledger 1/32=<path> 1/64=<path> \
        --out h_vs_ca.pdf [--csv h_vs_ca.csv]
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Computer Modern Roman"]
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "postProcess" / "reference"
RTUBE = 0.7


def aussillous_quere(ca):
    x = 1.34 * np.asarray(ca) ** (2.0 / 3.0)
    return x / (1.0 + 2.5 * x)


def bretherton(ca):
    return 1.34 * np.asarray(ca) ** (2.0 / 3.0)


def read_ledger(path: Path):
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = [r for r in data["results"] if r.get("status") == "SUCCESS"]
    ca_b = np.array([float(r["Ca_b"]) for r in rows])
    h = np.array([float(r["h_over_R0"]) for r in rows])
    order = np.argsort(ca_b)
    return ca_b[order], h[order], [rows[i] for i in order]


def read_csv(path: Path, ca_key: str, h_key: str, *, h_is_over_rtube: bool):
    with path.open(encoding="utf-8") as fp:
        rows = list(csv.DictReader(line for line in fp if not line.startswith("#")))
    ca = np.array([float(r[ca_key]) for r in rows])
    h = np.array([float(r[h_key]) for r in rows]) * (RTUBE if h_is_over_rtube else 1.0)
    order = np.argsort(ca)
    return ca[order], h[order]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ledger", nargs="+", required=True,
                   help="label=path pairs, e.g. 1/32=/run/root/continuation.json")
    p.add_argument("--out", type=Path, default=Path("h_vs_ca.pdf"))
    p.add_argument("--csv", type=Path, default=None, help="write the plotted rows")
    p.add_argument("--no-fem", action="store_true")
    args = p.parse_args(argv)

    fig, ax = plt.subplots(figsize=(12, 10))
    ca_law = np.logspace(-4, 1.2, 300)
    ax.plot(ca_law, bretherton(ca_law) * RTUBE, "--", color="0.55", lw=2.5,
            label=r"Bretherton, $1.34\,Ca_b^{2/3}$")
    ax.plot(ca_law, aussillous_quere(ca_law) * RTUBE, "-", color="0.25", lw=2.5,
            label=r"Aussillous--Qu\'er\'e")

    if not args.no_fem:
        ca, h = read_csv(REF / "fem-pr1-film-sweep.csv", "Ca_b", "H", h_is_over_rtube=True)
        ax.plot(ca, h, "s", ms=9, mfc="none", mec="tab:green", mew=2, label=r"FEM, void, $La=0$")
    ca, h = read_csv(REF / "fem-pr1-companion-basilisk.csv", "Ca_b", "H", h_is_over_rtube=True)
    ax.plot(ca, h, "^", ms=10, mfc="none", mec="tab:red", mew=2, label=r"lab frame, parent campaign")

    markers = ["o", "D", "P", "X"]
    colours = ["tab:blue", "tab:orange", "tab:purple", "tab:brown"]
    out_rows = []
    for k, spec in enumerate(args.ledger):
        label, _, path = spec.partition("=")
        ca_b, h, rows = read_ledger(Path(path))
        if ca_b.size == 0:
            continue
        ax.plot(ca_b, h, markers[k % 4], ms=11, color=colours[k % 4], mec="k", mew=1,
                zorder=4, label=rf"co-moving, $\Delta/R_0 = {label}$")
        for r in rows:
            out_rows.append({"delta": label, "Ca_in": r["Ca_in"], "Ca_b": r["Ca_b"],
                             "h_over_R0": r["h_over_R0"], "h_over_Rtube": r["h_over_Rtube"],
                             "minCells": r.get("minCells"), "route": r.get("route"),
                             "case_dir": r.get("case_dir")})

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$Ca_b = \mu U_b/\sigma$", fontsize=40, labelpad=15)
    ax.set_ylabel(r"$h/R_0$", fontsize=40, labelpad=15)
    ax.tick_params(which="both", direction="out", width=3, labelsize=30, pad=10)
    ax.tick_params(which="major", length=12)
    ax.tick_params(which="minor", length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(3)
    ax.minorticks_on()
    ax.legend(fontsize=22, frameon=False, loc="lower right")
    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, bbox_inches="tight", dpi=300, pad_inches=0.1)
    plt.close(fig)

    if args.csv:
        with args.csv.open("w", encoding="utf-8", newline="") as fp:
            w = csv.DictWriter(fp, fieldnames=list(out_rows[0].keys()) if out_rows else ["delta"])
            w.writeheader()
            w.writerows(out_rows)
    print(f"wrote {args.out} with {len(out_rows)} co-moving points")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
