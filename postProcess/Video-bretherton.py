#!/usr/bin/env python3
"""Render axisymmetric frames of a bretherton.c case and encode a video.

The tube is drawn about its axis with a different field in each half:
velocity magnitude above the axis, viscous dissipation below it. The
interface is mirrored so the bubble reads as one object.

Field data comes from ``postProcess/getData`` and the interface from
``postProcess/getFacets``; both are Basilisk helpers built from the same
snapshot, so the render never re-derives physics in Python.

Note on typography: this uses mathtext rather than ``text.usetex``.
LaTeX spawns a subprocess per unique string and deadlocks under
``multiprocessing``, which is exactly how frames are generated here.

Usage:
    python3 postProcess/Video-bretherton.py <case-dir> [--out video.mp4]
                                            [--cpus N] [--fps 20]
"""

import argparse
import os
import subprocess as sp
import sys
from functools import partial
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["mathtext.fontset"] = "cm"

HERE = os.path.dirname(os.path.abspath(__file__))
GETDATA = os.path.join(HERE, "getData")
GETFACETS = os.path.join(HERE, "getFacets")


def read_header(case_dir):
    """Return the run parameters recorded in the case log header."""
    log = [f for f in os.listdir(case_dir) if f.endswith("-log")]
    if not log:
        raise SystemExit(f"no case log in {case_dir}")
    with open(os.path.join(case_dir, log[0])) as fh:
        first = fh.readline()
    params = {}
    for item in first.lstrip("# ").strip().split(","):
        bits = item.split()
        if len(bits) == 2:
            params[bits[0]] = bits[1]
    return params


def facets(snapshot):
    """Interface segments as an (N, 2, 2) array of endpoint pairs."""
    out = sp.run([GETFACETS, snapshot], capture_output=True, text=True)
    segs, cur = [], []
    for line in out.stdout.splitlines():
        if not line.strip():
            if len(cur) == 2:
                segs.append(cur)
            cur = []
            continue
        a, b = line.split()
        cur.append((float(a), float(b)))
    if len(cur) == 2:
        segs.append(cur)
    return np.array(segs) if segs else np.empty((0, 2, 2))


def fields(snapshot, xmin, xmax, rmax, ny, muR):
    """Sample cs, f, log10 dissipation and |u| on a uniform grid."""
    out = sp.run([GETDATA, snapshot, str(xmin), "0", str(xmax), str(rmax),
                  str(ny), str(muR)], capture_output=True, text=True)
    rows = [r.split() for r in out.stderr.splitlines() if r.strip()]
    if not rows:
        return None
    d = np.array(rows, dtype=float)
    ncol = d.shape[1]          # getData's column count, not a fixed literal
    nyy = len(np.unique(d[:, 1]))
    nxx = len(d) // nyy
    d = d[: nxx * nyy].reshape(nxx, nyy, ncol)
    return d


def render(item, case, limits, outdir):
    idx, snapshot, tval = item
    dest = os.path.join(outdir, f"frame-{idx:05d}.png")
    if os.path.exists(dest):
        return dest

    seg = facets(snapshot)
    if len(seg) == 0:
        return None
    xlo = seg[:, :, 0].min() - 1.0
    xhi = seg[:, :, 0].max() + 1.0

    d = fields(snapshot, xlo, xhi, case["Rtube"], case["ny"], case["muR"])
    if d is None:
        return None
    X, Y = d[:, :, 0], d[:, :, 1]
    cs, diss, vel = d[:, :, 2], d[:, :, 4], d[:, :, 5]

    # Outside the fluid there is no field to show, only the wall.
    solid = ~(cs > 0.5)
    diss = np.where(solid, np.nan, diss)
    vel = np.where(solid, np.nan, vel)

    fig, ax = plt.subplots(figsize=(16, 5))
    fig.subplots_adjust(left=0.02, right=0.90, top=0.86, bottom=0.06)
    ext = [X.min(), X.max(), 0, case["Rtube"]]
    ax.imshow(vel.T, extent=ext, origin="lower", aspect="equal",
              cmap="Blues", vmin=limits["vmin"], vmax=limits["vmax"])
    ext_m = [X.min(), X.max(), -case["Rtube"], 0]
    ax.imshow(np.flipud(diss.T), extent=ext_m, origin="lower", aspect="equal",
              cmap="hot_r", vmin=limits["dmin"], vmax=limits["dmax"])

    lc = LineCollection(seg, colors="black", linewidths=2.0)
    ax.add_collection(lc)
    ax.add_collection(LineCollection(seg * np.array([1, -1]),
                                     colors="black", linewidths=2.0))

    for sgn in (1, -1):
        ax.plot([X.min(), X.max()], [sgn * case["Rtube"]] * 2,
                color="dimgrey", lw=4, solid_capstyle="butt")
    ax.axhline(0.0, color="grey", lw=0.8, ls=(0, (6, 6)))

    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(-case["Rtube"] * 1.05, case["Rtube"] * 1.05)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(rf"$Ca = {case['Ca']}$,  MAXlevel {case['MAXlevel']},  "
                 rf"$t = {tval:.2f}$", fontsize=22, pad=16)

    imv = plt.cm.ScalarMappable(cmap="Blues",
                                norm=plt.Normalize(limits["vmin"], limits["vmax"]))
    imd = plt.cm.ScalarMappable(cmap="hot_r",
                                norm=plt.Normalize(limits["dmin"], limits["dmax"]))
    cav = fig.add_axes([0.92, 0.53, 0.014, 0.34])
    cad = fig.add_axes([0.92, 0.13, 0.014, 0.34])
    cbv = fig.colorbar(imv, cax=cav)
    cbd = fig.colorbar(imd, cax=cad)
    cbv.set_label(r"$|u|$", fontsize=18, labelpad=10)
    cbd.set_label(r"$\log_{10}(\mu\,D\!:\!D)$", fontsize=18, labelpad=10)
    for cb in (cbv, cbd):
        cb.ax.tick_params(labelsize=13)

    fig.savefig(dest, dpi=110, pad_inches=0.08)
    plt.close(fig)
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--ny", type=int, default=140)
    args = ap.parse_args()

    case_dir = os.path.abspath(args.case_dir)
    p = read_header(case_dir)
    case = {"Ca": float(p["Ca"]), "Rtube": float(p["Rtube"]),
            "muR": float(p["muR"]), "MAXlevel": int(p["MAXlevel"]),
            "ny": args.ny}

    snapdir = os.path.join(case_dir, "intermediate")
    snaps = sorted(os.listdir(snapdir),
                   key=lambda s: float(s.split("-")[1]))
    items = [(i, os.path.join(snapdir, s), float(s.split("-")[1]))
             for i, s in enumerate(snaps)]
    if not items:
        raise SystemExit("no snapshots")

    # Fixed colour limits across the video, from a sample of frames, so
    # brightness changes mean physics rather than rescaling.
    probe = items[:: max(1, len(items) // 8)][:8]
    vv, dd = [], []
    for _, snap, _t in probe:
        seg = facets(snap)
        if len(seg) == 0:
            continue
        d = fields(snap, seg[:, :, 0].min() - 1.0, seg[:, :, 0].max() + 1.0,
                   case["Rtube"], case["ny"], case["muR"])
        if d is None:
            continue
        good = d[:, :, 2] > 0.5
        vv.append(np.nanpercentile(d[:, :, 5][good], 99.5))
        # getData writes -10 where the dissipation is identically zero.
        # Including that sentinel in the range would compress every real
        # value into the top of the colormap.
        dvals = d[:, :, 4][good]
        dvals = dvals[dvals > -9.9]
        if dvals.size:
            dd.append(np.nanpercentile(dvals, [2, 99.5]))
    if not vv:
        raise SystemExit("could not establish colour limits")
    if not dd:
        raise SystemExit("no finite dissipation values found")
    dmax = float(np.max([x[1] for x in dd]))
    dmin = float(np.min([x[0] for x in dd]))
    # Five decades is enough to show the film and the caps together; more
    # than that and the film detail washes out.
    dmin = max(dmin, dmax - 5.0)
    limits = {"vmin": 0.0, "vmax": float(np.max(vv)),
              "dmin": dmin, "dmax": dmax}
    print(f"colour limits: |u| 0 -> {limits['vmax']:.4g}, "
          f"log10 diss {limits['dmin']:.3g} -> {limits['dmax']:.3g}", flush=True)

    outdir = os.path.join(case_dir, "frames")
    os.makedirs(outdir, exist_ok=True)
    with Pool(args.cpus) as pool:
        made = pool.map(partial(render, case=case, limits=limits,
                                outdir=outdir), items)
    made = [m for m in made if m]
    print(f"rendered {len(made)} of {len(items)} frames", flush=True)

    out = args.out or os.path.join(case_dir,
                                   f"bretherton-Ca{case['Ca']}-ml{case['MAXlevel']}.mp4")
    cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
           "-pattern_type", "glob", "-i", os.path.join(outdir, "frame-*.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", out]
    r = sp.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-1500:], file=sys.stderr)
        raise SystemExit("ffmpeg failed")
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
