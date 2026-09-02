#!/usr/bin/env python3
"""Render axisymmetric frames of a bretherton.c case and encode a video.

Two rows: the whole domain on top, and a window that travels with the
bubble below. In both, the tube is drawn about its axis with a
different field in each half -- velocity magnitude above the axis,
viscous dissipation below it -- and the interface is mirrored so the
bubble reads as one object.

The travelling window has a **fixed width**. That is not cosmetic: the
embedded tube wall is rigid, so it must occupy the same pixels in every
frame. An earlier version padded a fixed margin around the interface,
so the window width tracked the bubble length; with ``aspect="equal"``
matplotlib then resized the axes box frame by frame and the wall
appeared to swell by 7% over a run. Fixed limits are what keep a rigid
wall looking rigid.

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


def read_ldomain(case_dir):
    """Domain length from the case parameter file."""
    pf = os.path.join(case_dir, "case.params")
    if os.path.exists(pf):
        for line in open(pf):
            line = line.split("#")[0].strip()
            if line.startswith("Ldomain"):
                return float(line.split("=")[1])
    return 16.0


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


def draw_panel(ax, d, seg, case, limits, xlo, xhi, equal):
    """One tube panel: |u| above the axis, dissipation below."""
    X = d[:, :, 0]
    cs, diss, vel = d[:, :, 2], d[:, :, 4], d[:, :, 5]
    solid = ~(cs > 0.5)
    diss = np.where(solid, np.nan, diss)
    vel = np.where(solid, np.nan, vel)

    Rt = case["Rtube"]
    ax.imshow(vel.T, extent=[X.min(), X.max(), 0, Rt], origin="lower",
              aspect="auto", cmap="Blues",
              vmin=limits["vmin"], vmax=limits["vmax"])
    ax.imshow(np.flipud(diss.T), extent=[X.min(), X.max(), -Rt, 0],
              origin="lower", aspect="auto", cmap="hot_r",
              vmin=limits["dmin"], vmax=limits["dmax"])

    if len(seg):
        ax.add_collection(LineCollection(seg, colors="black", linewidths=1.8))
        ax.add_collection(LineCollection(seg * np.array([1, -1]),
                                         colors="black", linewidths=1.8))

    # The wall is rigid: draw it at +/-Rtube, never at a data extreme.
    for sgn in (1, -1):
        ax.plot([xlo, xhi], [sgn * Rt] * 2, color="dimgrey", lw=5,
                solid_capstyle="butt", zorder=5)
    ax.axhline(0.0, color="grey", lw=0.7, ls=(0, (6, 6)), zorder=4)

    # Fixed limits every frame -> fixed pixels for the wall.
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(-Rt * 1.06, Rt * 1.06)
    if equal:
        # adjustable="box" keeps the limits I set and shrinks the axes
        # box instead. Because the window width is fixed, that shrink is
        # identical in every frame, so the wall holds its pixels.
        # adjustable="datalim" would silently rewrite the y limits.
        ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def render(item, case, limits, outdir):
    idx, snapshot, tval = item
    dest = os.path.join(outdir, f"frame-{idx:05d}.png")
    if os.path.exists(dest):
        return dest

    seg = facets(snapshot)
    if len(seg) == 0:
        return None

    Rt, Ld, W = case["Rtube"], case["Ldomain"], case["window"]

    # Travelling window: fixed width, centred on the bubble, clamped to
    # the domain so the width can never change at the ends.
    centre = 0.5 * (seg[:, :, 0].min() + seg[:, :, 0].max())
    xlo = min(max(centre - W / 2.0, 0.0), Ld - W)
    xhi = xlo + W

    d_full = fields(snapshot, 0.0, Ld, Rt, case["ny_full"], case["muR"])
    d_win = fields(snapshot, xlo, xhi, Rt, case["ny"], case["muR"])
    if d_full is None or d_win is None:
        return None

    # The lower panel holds true aspect, so its height is set by the
    # window width; the figure is sized to that rather than leaving the
    # shrunk axes floating in whitespace.
    panel_w = 17.0 * (0.88 - 0.03)
    win_h = panel_w * (2.0 * Rt * 1.06) / W
    fig_h = win_h + 1.55 + 1.75
    fig = plt.figure(figsize=(17, fig_h))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.55, win_h],
                          left=0.03, right=0.88,
                          top=1.0 - 1.00 / fig_h, bottom=0.30 / fig_h,
                          hspace=0.55 / max(win_h, 0.5))
    ax_full = fig.add_subplot(gs[0])
    ax_win = fig.add_subplot(gs[1])

    draw_panel(ax_full, d_full, seg, case, limits, 0.0, Ld, equal=False)
    draw_panel(ax_win, d_win, seg, case, limits, xlo, xhi, equal=True)

    # Mark where the travelling window sits within the whole domain.
    for xv in (xlo, xhi):
        ax_full.plot([xv, xv], [-Rt * 1.06, Rt * 1.06], color="black",
                     lw=1.2, ls="--", zorder=6)

    ax_full.set_title(rf"$L = {Ld:g}$", fontsize=14, pad=5)
    ax_win.set_title(rf"window ${W:g}R$", fontsize=14, pad=5)
    fig.suptitle(rf"$Ca = {case['Ca']}$,  MAXlevel {case['MAXlevel']},  "
                 rf"$t = {tval:.2f}$", fontsize=21,
                 y=1.0 - 0.20 / fig_h)

    imv = plt.cm.ScalarMappable(cmap="Blues",
                                norm=plt.Normalize(limits["vmin"], limits["vmax"]))
    imd = plt.cm.ScalarMappable(cmap="hot_r",
                                norm=plt.Normalize(limits["dmin"], limits["dmax"]))
    cbv = fig.colorbar(imv, cax=fig.add_axes([0.90, 0.54, 0.012, 0.30]))
    cbd = fig.colorbar(imd, cax=fig.add_axes([0.90, 0.14, 0.012, 0.30]))
    cbv.set_label(r"$|u|$", fontsize=17, labelpad=9)
    cbd.set_label(r"$\log_{10}(\mu\,D\!:\!D)$", fontsize=17, labelpad=9)
    for cb in (cbv, cbd):
        cb.ax.tick_params(labelsize=12)

    fig.savefig(dest, dpi=105)
    plt.close(fig)
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case_dir")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--ny", type=int, default=140,
                    help="radial samples in the travelling window")
    ap.add_argument("--ny-full", type=int, default=90,
                    help="radial samples across the full domain")
    ap.add_argument("--window", type=float, default=8.0,
                    help="travelling window width in units of R (fixed)")
    args = ap.parse_args()

    case_dir = os.path.abspath(args.case_dir)
    p = read_header(case_dir)
    Ld = read_ldomain(case_dir)
    case = {"Ca": float(p["Ca"]), "Rtube": float(p["Rtube"]),
            "muR": float(p["muR"]), "MAXlevel": int(p["MAXlevel"]),
            "ny": args.ny, "ny_full": args.ny_full,
            "Ldomain": Ld, "window": min(args.window, Ld)}

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
        c = 0.5 * (seg[:, :, 0].min() + seg[:, :, 0].max())
        wlo = min(max(c - case["window"] / 2.0, 0.0),
                  case["Ldomain"] - case["window"])
        d = fields(snap, wlo, wlo + case["window"],
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
