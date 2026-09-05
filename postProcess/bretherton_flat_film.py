#!/usr/bin/env python3
"""Flat-film thickness from saved snapshots, without re-running anything.

`bretherton.c` logs ``bFilm = Rtube - max_y(interface)``: the *thinnest*
point anywhere on the bubble. That minimum sits in the rear-meniscus
ripple, not in the uniform film, so it understates the flat film by 20
to 45% and increasingly so as Ca falls. Bretherton's and Aussillous and
Quere's laws describe the *flat* film, so the comparison needs that.

Both quantities are recoverable from the interface geometry alone, so
this reads the ``intermediate/snapshot-*`` series through
``postProcess/getFacets`` and reports:

- ``b_global`` = Rtube - max(r), reproducing the solver's ``bFilm``;
- ``b_flat``   = Rtube - median(r) over the central part of the bubble
  body, trimming a fraction of the length at each meniscus.

The trim fraction is the one free choice. Vary it with ``--trim`` to
check it does not matter: between 0.20 and 0.40 the answer moves by less
than 0.5%, well under the effect being measured.

Usage:
    python3 postProcess/bretherton_flat_film.py <case-dir> [...]
                [--trim 0.30] [--window 0.25] [--cpus 4] [--out flat.csv]
"""

import argparse
import glob
import os
import subprocess as sp
from functools import partial
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GETFACETS = os.path.join(HERE, "getFacets")


def header(case_dir):
    logs = [f for f in os.listdir(case_dir) if f.endswith("-log")]
    if not logs:
        raise SystemExit(f"no case log in {case_dir}")
    with open(os.path.join(case_dir, logs[0])) as fh:
        first = fh.readline()
    out = {}
    for item in first.lstrip("# ").strip().split(","):
        bits = item.split()
        if len(bits) == 2:
            out[bits[0]] = bits[1]
    return out


def measure(snapshot, Rtube, trim):
    """Return (t, xTipR, xTipF, b_global, b_flat) for one snapshot."""
    t = float(os.path.basename(snapshot).split("-")[-1])
    res = sp.run([GETFACETS, snapshot], capture_output=True, text=True)
    pts = [l.split() for l in res.stdout.splitlines() if len(l.split()) == 2]
    if len(pts) < 8:
        return None
    a = np.array(pts, dtype=float)
    x, r = a[:, 0], a[:, 1]
    xr, xf = x.min(), x.max()
    L = xf - xr
    if L <= 0:
        return None
    mid = (x > xr + trim * L) & (x < xf - trim * L)
    if mid.sum() < 4:
        return None
    return (t, xr, xf, Rtube - r.max(), Rtube - float(np.median(r[mid])))


def analyse(case_dir, trim, window, cpus):
    p = header(case_dir)
    Rtube, Ca = float(p["Rtube"]), float(p["Ca"])
    snaps = sorted(glob.glob(os.path.join(case_dir, "intermediate", "snapshot-*")),
                   key=lambda s: float(s.split("-")[-1]))
    if not snaps:
        raise SystemExit(f"no snapshots in {case_dir}")
    with Pool(cpus) as pool:
        rows = pool.map(partial(measure, Rtube=Rtube, trim=trim), snaps)
    rows = [r for r in rows if r]
    if len(rows) < 4:
        raise SystemExit(f"too few usable snapshots in {case_dir}")
    d = np.array(rows)
    t, xr, xf, bg, bf = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]

    # Tip capillary number from the front-tip trajectory over the final
    # window, matching how bretherton_film.py defines Ca_b.
    n = max(3, int(len(t) * window))
    Ca_b = float(np.polyfit(t[-n:], xf[-n:], 1)[0])
    b_flat = float(np.mean(bf[-n:]))
    b_glob = float(np.mean(bg[-n:]))

    bret = 1.34 * Ca_b ** (2 / 3)
    aq = bret / (1 + 3.35 * Ca_b ** (2 / 3))
    return {"case": os.path.basename(os.path.normpath(case_dir)), "Ca": Ca,
            "Ca_b": Ca_b, "Rtube": Rtube,
            "b_flat_over_R": b_flat / Rtube, "b_global_over_R": b_glob / Rtube,
            "ratio": b_flat / b_glob if b_glob else float("nan"),
            "bretherton": bret, "aq": aq,
            "dev_flat_aq": 100 * (b_flat / Rtube - aq) / aq,
            "dev_global_aq": 100 * (b_glob / Rtube - aq) / aq,
            "series": d}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dirs", nargs="+")
    ap.add_argument("--trim", type=float, default=0.30,
                    help="fraction of bubble length trimmed at each meniscus")
    ap.add_argument("--window", type=float, default=0.25,
                    help="trailing fraction of snapshots averaged")
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--out", default=None, help="write a summary CSV")
    args = ap.parse_args()
    if not 0.0 < args.trim < 0.5:
        ap.error("--trim must lie in (0, 0.5)")
    if not 0.0 < args.window <= 1.0:
        ap.error("--window must lie in (0, 1]")

    results = [analyse(c, args.trim, args.window, args.cpus)
               for c in args.case_dirs]
    print(f"{'case':>6} {'Ca_b':>9} {'b_flat/R':>9} {'b_glob/R':>9} "
          f"{'flat/glob':>9} {'AQ':>9} {'dev flat':>9} {'dev glob':>9}")
    for r in results:
        print(f"{r['case']:>6} {r['Ca_b']:>9.5f} {r['b_flat_over_R']:>9.5f} "
              f"{r['b_global_over_R']:>9.5f} {r['ratio']:>9.3f} {r['aq']:>9.5f} "
              f"{r['dev_flat_aq']:>8.2f}% {r['dev_global_aq']:>8.2f}%")
    if args.out:
        with open(args.out, "w") as fh:
            fh.write("# flat-film thickness recovered from saved snapshots\n")
            fh.write(f"# trim={args.trim} window={args.window}\n")
            fh.write("case,Ca,Ca_b,b_flat_over_R,b_global_over_R,ratio,"
                     "bretherton,aq,dev_flat_aq,dev_global_aq\n")
            for r in results:
                fh.write(f"{r['case']},{r['Ca']},{r['Ca_b']:.8g},"
                         f"{r['b_flat_over_R']:.8g},{r['b_global_over_R']:.8g},"
                         f"{r['ratio']:.6g},{r['bretherton']:.8g},{r['aq']:.8g},"
                         f"{r['dev_flat_aq']:.4g},{r['dev_global_aq']:.4g}\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
