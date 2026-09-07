"""Render the approved flat-film comparison for the Newtonian report.

Input: ../data/baseline-film.csv. Outputs: PDF and PNG beside this script.
The Taylor reference is the correlation in Aussillous & Quere (2000).
"""

from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
import numpy as np


ROOT = Path(__file__).resolve().parent
with (ROOT.parent / "data/baseline-film.csv").open() as source:
    rows = list(csv.DictReader(line for line in source if not line.startswith("#")))

ca = np.array([float(row["Ca_b"]) for row in rows])
film = np.array([float(row["b_flat_over_R"]) for row in rows])
assert len(ca) == 4 and np.all(ca > 0) and np.all(film > 0)


def bretherton(x):
    return 1.34 * np.asarray(x) ** (2 / 3)


def taylor(x):
    b = bretherton(x)
    return b / (1 + 2.5 * b)


b_ref, t_ref = bretherton(ca), taylor(ca)
assert np.allclose(t_ref, [float(row["taylor_aq_b_over_R"]) for row in rows], rtol=1e-6)
curve_ca = np.geomspace(1e-4, 1, 501)

matplotlib.rcParams.update({
    "font.family": "serif", "font.serif": ["Computer Modern Roman"],
    "text.usetex": True, "text.latex.preamble": r"\usepackage{amsmath}",
})
fig, (ax, ratio_ax) = plt.subplots(1, 2, figsize=(17, 8),
                                 gridspec_kw={"width_ratios": [1.35, 1]})
orange, blue = "#D55E00", "#0072B2"
ax.plot(curve_ca, taylor(curve_ca), color=orange, lw=3,
        label=r"Taylor law (Aussillous--Qu\'er\'e fit)")
# This is a limiting formula, not a high-Ca model; do not extend it to Ca=1.
small_ca = np.geomspace(1e-4, 0.1, 350)
ax.plot(small_ca, bretherton(small_ca), "--", color="0.4", lw=2.7,
        label=r"Bretherton: $1.34\,Ca_b^{2/3}$")
ax.scatter(ca, film, s=170, color=blue, edgecolor="white", linewidth=1.2,
           zorder=5, label="Bubble simulations (flat film)")
ax.set(xscale="log", yscale="log", xlim=(8e-5, 1.2), ylim=(.0018, .5))
ax.set_ylabel(r"$b_{\mathrm{flat}}/R_{\mathrm{tube}}$", fontsize=32, labelpad=12)
ax.legend(loc="upper left", fontsize=19, frameon=False, handlelength=2.2)
ax.text(.03, .05, r"$La_v=1,\quad \mu_g/\mu_c=0.01,\quad \rho_g/\rho_c=10^{-3}$",
        transform=ax.transAxes, fontsize=18)

ratio_ax.axhline(1, color=orange, lw=2.5, label="Taylor correlation")
ratio_ax.scatter(ca, film/t_ref, s=170, color=blue, edgecolor="white",
                 linewidth=1.2, zorder=5)
ratio_ax.set(xscale="log", xlim=(.004, .085), ylim=(.975, 1.009))
ratio_ax.set_ylabel(r"$b_{\mathrm{flat}}/b_{\mathrm{Taylor}}$", fontsize=30, labelpad=12)
ratio_ax.xaxis.set_major_locator(FixedLocator([.005, .01, .02, .05]))
ratio_ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
ratio_ax.xaxis.set_minor_formatter(NullFormatter())
ratio_ax.yaxis.set_major_locator(FixedLocator([.98, .99, 1.00]))
ratio_ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.2f}"))
ratio_ax.text(.06, .08, "Completed cases only", transform=ratio_ax.transAxes,
              fontsize=19)
for x, y in zip(ca, film/t_ref):
    ratio_ax.annotate(f"${100*(y-1):+.2f}\\%$", (x, y),
                      xytext=(0, 13), textcoords="offset points",
                      ha="center", fontsize=17)

for index, axes in enumerate((ax, ratio_ax)):
    axes.set_xlabel(r"$Ca_b=\mu_c U_b/\sigma$", fontsize=31, labelpad=12)
    axes.tick_params(which="both", direction="out", width=2, labelsize=23, pad=7)
    axes.tick_params(which="major", length=9)
    axes.tick_params(which="minor", length=4.5)
    axes.text(0, 1.03, f"({chr(97+index)})", transform=axes.transAxes, fontsize=24)
    for spine in axes.spines.values():
        spine.set_linewidth(2)
fig.tight_layout(w_pad=3)
fig.savefig(ROOT / "film-thickness-vs-ca.pdf", dpi=300,
            bbox_inches="tight", pad_inches=.15)
fig.savefig(ROOT / "film-thickness-vs-ca.png", dpi=300,
            bbox_inches="tight", pad_inches=.15)
plt.close(fig)
print("Rendered four measured points, Taylor correlation, and Bretherton asymptote.")
