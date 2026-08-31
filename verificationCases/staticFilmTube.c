/**
# staticFilmTube.c

**Verification** of curvature and property-jump errors for a thin
interface--wall separation inside an embedded tube (evidence source:
exact Young--Laplace solution of the implemented equations).

## Claim

At production property ratios, the discrete solver reproduces the exact
static solution of the equations it implements -- a bubble at rest with
$\Delta p = 2\sigma/R_d$ and $\mathbf{u} \equiv 0$ -- to within a stated
tolerance, and the error grows as the interface is brought to within a
few cells of the embedded wall.

## Comparator

The comparator is an *exact solution of the implemented equations*, not
an experiment or a correlation: a sphere of radius $R_d$ held in
equilibrium by surface tension in a quiescent, gravity-free,
axisymmetric domain. For that state the momentum equation is satisfied
by $\mathbf{u} = 0$ with a piecewise-uniform pressure, so the exact
solution is independent of how close the sphere sits to the tube wall
and of the wall's own no-slip condition. This is therefore verification:
every non-zero velocity is a spurious (parasitic) current, and every
deviation of the measured jump from $2\sigma/R_d$ is discretisation
error.

The third coupling documented in `src-local/embed-vof-tube.h` is the one
under test: `heights.h`/`curvature.h` have no knowledge of `cs`, so the
film between the interface and the embedded wall must remain resolved by
several cells. This case measures the price of violating that.

## Motivation

The production campaign (`simulationCases/bretherton.c`) runs a bubble
at $R_{tube} = 0.7$, $L_{domain} = 16$, `MAXlevel = 12`, i.e.
$\Delta = 3.906\times 10^{-3}$. The Bretherton estimate
$b = 1.34\,Ca^{2/3} R_{tube}$ gives roughly 32 cells of film at
$Ca = 0.05$ but only about 3.8 cells at $Ca = 0.002$, below the
repository's own four-cell warning threshold. No existing verification
case exercises curvature evaluation at a few-cell interface--wall gap,
and none carries any density or viscosity contrast: `laplaceEmbedTube.c`
uses matched properties and so cannot probe the property-jump path at
all. This case closes both gaps.

Unlike `laplaceEmbedTube.c` this case does **not** define `FILTERED`, so
the density and viscosity fields are constructed exactly as in
`simulationCases/bretherton.c`.

## Geometry and sampling regions

Axisymmetric domain $L_0 = 2$, uniform grid $N = 256$, hence
$\Delta = 7.8125\times 10^{-3}$. The tube wall is embedded at
$y = R_{tube} = 0.7$ (fluid below). A spherical bubble of radius $R_d$
is centred on the axis at $x = X_d = 0.75$; it never touches the wall,
so the equilibrium shape is exactly a sphere and the interface--wall gap
is $R_{tube} - R_d$.

Because $R_d \approx 0.65$ leaves essentially no fluid at
$r > 1.5 R_d$ inside the tube, the continuous phase cannot be sampled
radially as in `laplaceEmbedTube.c`. It is sampled **axially displaced**
instead:

- *inside* (dispersed phase): full cells (`cs >= 1`) with
  $r = \sqrt{(x-X_d)^2 + y^2} < R_d/2$;
- *outside* (continuous phase): full cells with $x > 1.6$, i.e. at least
  23 cells beyond the bubble's far tip at $x = X_d + R_d \approx 1.42$.

Both regions are counted, and their phase purity is asserted
($f > 1 - 10^{-10}$ inside, $f < 10^{-10}$ outside); an empty, small or
impure region is a hard failure, since it would silently corrupt the
measured jump.

## Gap sweep

$R_d$ is chosen to give gaps of exactly 4, 6 and 8 cells:
$R_d = R_{tube} - n\Delta$ for $n = 4, 6, 8$, i.e.
$R_d = 0.66875,\ 0.653125,\ 0.6375$.

## Tolerance control

The run to $t_{end} = 2$ is roughly three visco-capillary times
$\mu_c R_d/\sigma \approx 0.65$, long enough for the parasitic-current
amplitude to be more than a transient. The 4-cell gap -- the headline
admissibility case -- is repeated with a ten-fold tighter Poisson
`TOLERANCE` ($10^{-6}$ against the baseline $10^{-5}$). If $\Delta p$
and $\max|\mathbf{u}|$ are unmoved, the measured error is discretisation
error and can be attributed to curvature and the property jump rather
than to an under-converged projection. "Unmoved" is defined
non-arbitrarily: for $\Delta p$, a relative shift below a tenth of the
gate tolerance; for $\max|\mathbf{u}|$, a relative shift below twice the
within-run coefficient of variation of that same statistic over
$t \ge t_{end}/2$, which is measured and reported rather than assumed.

A third point is added at the tolerance the production campaign actually
uses, `TOLERANCE` $= 10^{-4}$, which is *looser* than the baseline here.
The three-point ladder $10^{-4}, 10^{-5}, 10^{-6}$ at the 4-cell gap
therefore says directly whether the projection tolerance or the film
resolution is the binding constraint on the production configuration.

The 6- and 8-cell gaps are not tolerance-controlled; that remains an
untested assumption.

## Pass predicate

The gate is applied at the **largest (8-cell) gap**, where the film is
comfortably resolved:

- relative $\Delta p$ error below 2%;
- sustained spurious-current amplitude below
  $2\times 10^{-4}\,\sigma/\mu_c$. This threshold is not arbitrary: the
  smallest production capillary number is $Ca = 0.002$, so the imposed
  velocity there is $0.002\,\sigma/\mu_c$ in code units, and the
  threshold demands that spurious currents stay below 10% of it. The
  gated statistic is the *peak* of $\max|\mathbf{u}|$ over the logged
  samples in the second half of the run, $t \ge t_{end}/2$, not the
  instantaneous value at $t_{end}$: the parasitic-current amplitude
  fluctuates by tens of percent between samples, so a snapshot is not a
  reproducible bound whereas a second-half peak is a conservative one;
- every sampling region non-empty, adequately populated and phase-pure.

Gating at 8 cells makes the 4-cell result *diagnostic* rather than
self-referential: if the solver is certified at a comfortably resolved
separation, a 4-cell failure is attributable to the film resolution and
not to a broken solver. The full table is printed unconditionally and is
the scientific product; gaps that miss the thresholds are flagged with
explicit advisory lines and, for the 4-cell gap, an admissibility
verdict. Machine-readable results go to `staticFilmTube-gap.csv`.

## Limits

This is a **static** bound. It constrains curvature and property-jump
error at a given interface--wall separation on a uniform grid, with a
spherical interface and a quiescent exact solution. It says nothing
directly about whether the *dynamic* Bretherton film at the same
separation is adequately resolved: the dynamic film is not spherical,
carries a lubrication pressure gradient, is advected, and in production
sits on an adapted tree rather than a uniform grid. A pass here is
necessary, not sufficient. Conversely, a failure here is a genuine
obstruction, since the dynamic problem is strictly harder.

Three further limits are worth naming explicitly.

First, whenever the tolerance ladder shows $\max|\mathbf{u}|$ moving
with `TOLERANCE`, the tabulated spurious-current amplitude is an *upper
bound set by the projection*, not a measurement of the curvature-driven
current. Under those conditions the velocity column bounds the solver
configuration and not the surface-tension discretisation, and the
$\Delta p$ column is the only quantity in the table that speaks to
curvature. The case prints which of the two situations holds rather than
leaving it to the reader.

Second, the $\Delta p$ errors reported here are comparable in magnitude
to the temporal fluctuation of $\Delta p$ within a single run, so the
*ordering* of the three gaps is not resolved; only their common
magnitude, and its change with $\Delta$, are meaningful.

Third, the gap sweep varies the interface--wall separation at fixed
$\Delta$, which necessarily also varies $R_d/\Delta$. The two effects are
not separated within a single grid; they are distinguished only by
comparing runs at different $N$.
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "embed-vof-tube.h"

/**
## Parameters

Code units $\sigma = \mu_c = R = 1$, matching `simulationCases/bretherton.c`.
*/

#define Rtube    0.7
#define Ldom     2.0
#define Ngrid    256
#define hgrid    (Ldom/(double)Ngrid)
#define Xd       0.75

#define tend     2.0
#define tlog     0.05
/** History buffer; must exceed `tend/tlog` (an integer constant
expression is required at file scope). */
#define NHIST    64

#define TOL_BASE 1e-5
#define TOL_TIGHT 1e-6
/** The projection tolerance actually used by
`simulationCases/bretherton.c`. */
#define TOL_PROD 1e-4
#define DTMAX    1e-2

/** Sampling regions. */
#define RIN_FRAC 0.5
#define XOUT     1.6
#define NMIN_IN  100
#define NMIN_OUT 100
#define FPURE    1e-10

/** A shift in $\Delta p$ under the tolerance control is material if it
is a tenth of the gate tolerance or more; below that the projection
tolerance cannot plausibly account for the measured curvature error. */
#define TOL_DP_TOLCTL 2e-3

/** Gate, applied at the 8-cell gap only (see the pass predicate above). */
#define GAP_GATE 8.
#define TOL_DP   2e-2
#define TOL_U    2e-4

/** Plateau criterion: the mean of $\max|u|$ over the last quarter of the
run no more than 5% above its mean over the second quarter. Windowed
means are used because the instantaneous amplitude fluctuates by tens of
percent between logged samples. */
#define PLATEAU_RATIO 1.05

#define NRUN 5

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

/**
## Run bookkeeping

`Rd` and `TOLERANCE` are reset before each call to `run()`; the results
of every run are collected here and reported together at the end.
*/

typedef struct {
  double gapcells, Rd, gap, tol;
  double dp, dpex, edp;
  double umax, umax_cut;
  double upeak, tpeak, usus, uwin1, uwin2, uratio, ucv;
  double Rdeff, dpeff;
  int nin, nout, nsteps;
  double fmin_in, fmax_out;
  bool pure, populated;
} Result;

static double Rd = 0.;
static int run_id = 0;
static Result res[NRUN];
static double thist[NHIST], uhist[NHIST];
static int nhist = 0;

int main()
{
  L0 = Ldom;
  N  = Ngrid;

  /**
  Fluid 1 (`f = 1`) is the dispersed bubble; fluid 2 (`f = 0`) is the
  continuous wetting phase. The production bubble defaults are
  `La = 1`, `muR = 0.01`, `rhoR = 0.001`, and the property construction
  below is copied verbatim from `simulationCases/bretherton.c`. */

  double La = 1., muR = 0.01, rhoR = 0.001;
  rho1 = La*rhoR; mu1 = muR;
  rho2 = La;      mu2 = 1.;
  f.sigma = 1.;

  DT = DTMAX;

  /** Three gaps at the baseline tolerance, then the 4-cell gap repeated
  at a ten-fold tighter and at the production projection tolerance. */

  double gapcells[NRUN] = {4., 6., 8., 4., 4.};
  double tols[NRUN]     = {TOL_BASE, TOL_BASE, TOL_BASE, TOL_TIGHT, TOL_PROD};

  for (run_id = 0; run_id < NRUN; run_id++) {
    Rd = Rtube - gapcells[run_id]*hgrid;
    TOLERANCE = tols[run_id];
    res[run_id].gapcells = gapcells[run_id];
    res[run_id].tol = TOLERANCE;
    nhist = 0;
    fprintf (stderr, "# run %d: gap = %g cells, Rd = %.6f, TOLERANCE = %g\n",
             run_id, gapcells[run_id], Rd, TOLERANCE);
    run();
  }

  /**
  ## Report
  */

  double delta = hgrid;
  printf ("== staticFilmTube: static bubble near an embedded tube wall ==\n");
  printf ("grid        : L0 = %g, N = %d, Delta = %.6e, Rtube = %g, Xd = %g\n",
          Ldom, Ngrid, delta, Rtube, Xd);
  printf ("properties  : rho1 = %g, mu1 = %g (dispersed) | "
          "rho2 = %g, mu2 = %g (continuous) | sigma = %g\n",
          rho1, mu1, rho2, mu2, f.sigma);
  printf ("run control : DT = %g, tend = %g, baseline TOLERANCE = %g, "
          "tight TOLERANCE = %g\n", DTMAX, tend, TOL_BASE, TOL_TIGHT);
  printf ("sampling    : inside r < %g*Rd about (%g,0); "
          "outside x > %g with cs >= 1\n\n", RIN_FRAC, Xd, XOUT);

  printf ("%-5s %-9s %-10s %-8s %-12s %-12s %-10s %-11s %-11s %-9s\n",
          "gap", "Rd", "TOLERANCE", "steps", "dp", "dp_exact", "rel.err",
          "max|u|", "max|u|_cut", "Rd_eff");
  printf ("%-5s %-9s %-10s %-8s %-12s %-12s %-10s %-11s %-11s %-9s\n",
          "[cell]", "[R]", "", "", "[sig/R]", "[sig/R]", "[-]",
          "[sig/mu]", "[sig/mu]", "[R]");
  for (int i = 0; i < NRUN; i++) {
    Result * r = &res[i];
    printf ("%-5.1f %-9.6f %-10.1e %-8d %-12.6e %-12.6e %-10.3e "
            "%-11.4e %-11.4e %-9.6f\n",
            r->gapcells, r->Rd, r->tol, r->nsteps, r->dp, r->dpex, r->edp,
            r->umax, r->umax_cut, r->Rdeff);
  }

  printf ("\nSampling regions and phase purity\n");
  printf ("%-5s %-10s %-8s %-9s %-14s %-14s\n",
          "gap", "TOLERANCE", "n_in", "n_out", "min(f) inside",
          "max(f) outside");
  for (int i = 0; i < NRUN; i++) {
    Result * r = &res[i];
    printf ("%-5.1f %-10.1e %-8d %-9d %-14.6e %-14.6e\n",
            r->gapcells, r->tol, r->nin, r->nout, r->fmin_in, r->fmax_out);
  }

  printf ("\nParasitic-current time history (has max|u| plateaued by tend?)\n");
  printf ("%-5s %-10s %-12s %-13s %-12s %-9s %-11s %-11s %-7s %-7s %s\n",
          "gap", "TOLERANCE", "u(tend)", "sustained", "peak", "t_peak",
          "mean 2nd qu", "mean 4th qu", "ratio", "cv", "verdict");
  for (int i = 0; i < NRUN; i++) {
    Result * r = &res[i];
    printf ("%-5.1f %-10.1e %-12.4e %-13.4e %-12.4e %-9.3f %-11.4e %-11.4e "
            "%-7.3f %-7.3f %s\n",
            r->gapcells, r->tol, r->umax, r->usus, r->upeak, r->tpeak,
            r->uwin1, r->uwin2, r->uratio, r->ucv,
            r->uratio <= PLATEAU_RATIO ?
            "plateaued/decaying" : "STILL GROWING at tend");
  }
  printf ("  \"sustained\" = peak of max|u| over logged samples with "
          "t >= tend/2; this is the gated statistic.\n");

  /**
  ### Tolerance control

  Runs 0 and 3 are the same 4-cell configuration at $10^{-5}$ and
  $10^{-6}$. A material shift would mean the reported error is
  contaminated by the projection tolerance and cannot be attributed to
  curvature.
  */

  double ddp = fabs (res[3].dp - res[0].dp)/fabs (res[0].dp);
  double du  = (res[0].usus > 0. ?
                fabs (res[3].usus - res[0].usus)/res[0].usus : HUGE);
  double ucv = (res[0].ucv > res[3].ucv ? res[0].ucv : res[3].ucv);
  bool dp_moved = (ddp >= TOL_DP_TOLCTL);
  bool u_moved  = (du >= 2.*ucv);

  printf ("\nTolerance ladder at the 4-cell gap "
          "(production TOLERANCE is %g)\n", TOL_PROD);
  printf ("  %-12s %-14s %-11s %-15s %-7s\n",
          "TOLERANCE", "dp", "rel.err", "sustained max|u|", "cv");
  int ladder[3] = {4, 0, 3};
  for (int k = 0; k < 3; k++) {
    Result * r = &res[ladder[k]];
    printf ("  %-12.1e %-14.6e %-11.3e %-15.4e %-7.3f\n",
            r->tol, r->dp, r->edp, r->usus, r->ucv);
  }

  printf ("\nTolerance control (4-cell gap, TOLERANCE %g -> %g)\n",
          TOL_BASE, TOL_TIGHT);
  printf ("  dp        : %.6e -> %.6e  (relative change %.3e; "
          "material if >= %.1e)\n",
          res[0].dp, res[3].dp, ddp, TOL_DP_TOLCTL);
  printf ("  max|u|    : %.6e -> %.6e  (sustained, relative change %.3e)\n",
          res[0].usus, res[3].usus, du);
  printf ("  scatter   : within-run coefficient of variation of max|u| over "
          "t >= tend/2 is %.3f (baseline) and %.3f (tight); a change is "
          "called material only above twice the larger, i.e. %.3f\n",
          res[0].ucv, res[3].ucv, 2.*ucv);
  printf ("  => dp %s; max|u| %s\n",
          dp_moved ? "MOVED MATERIALLY" : "unmoved",
          u_moved ? "MOVED MATERIALLY" : "unmoved (within its own scatter)");
  printf ("  => the measured pressure jump %s\n", dp_moved ?
          "cannot be attributed to curvature: it is contaminated by the "
          "projection tolerance"
          : "is discretisation error, not projection-tolerance error, and "
          "may be attributed to curvature and the property jump");
  printf ("  => the spurious-current amplitude %s\n", u_moved ?
          "is NOT a curvature bound: it is set by the projection tolerance, "
          "so the tabulated max|u| is an upper bound on the curvature-driven "
          "current and not a measurement of it"
          : "is insensitive to the projection tolerance and may be read as a "
          "curvature-driven amplitude");

  /**
  ### Advisories and gate
  */

  printf ("\nAdvisories\n");
  int nadv = 0;
  for (int i = 0; i < NRUN; i++) {
    Result * r = &res[i];
    if (r->edp >= TOL_DP) {
      printf ("  ADVISORY: gap %.1f cells (TOLERANCE %g): relative dp error "
              "%.3e exceeds %.3e\n", r->gapcells, r->tol, r->edp, TOL_DP);
      nadv++;
    }
    if (r->usus >= TOL_U) {
      printf ("  ADVISORY: gap %.1f cells (TOLERANCE %g): sustained max|u| "
              "%.4e exceeds %.4e sigma/mu_c\n",
              r->gapcells, r->tol, r->usus, TOL_U);
      nadv++;
    }
    if (r->uratio > PLATEAU_RATIO) {
      printf ("  ADVISORY: gap %.1f cells (TOLERANCE %g): max|u| has not "
              "plateaued by tend (mean over the last quarter is %.3f times "
              "the mean over the second quarter); the reported amplitude is a "
              "lower bound, not a converged value\n",
              r->gapcells, r->tol, r->uratio);
      nadv++;
    }
  }
  if (!nadv)
    printf ("  none\n");

  /** The 4-cell gap answers the admissibility question for the
  $Ca = 0.002$ production case at `MAXlevel = 12`, where the Bretherton
  film is about 3.8 cells thick. It is reported here, prominently, and is
  deliberately *not* used to relax anything. */

  bool ok4 = (res[0].edp < TOL_DP && res[0].usus < TOL_U);
  printf ("\nADMISSIBILITY (4-cell interface--wall gap, the Ca = 0.002 "
          "production condition at MAXlevel 12)\n");
  printf ("  relative dp error %.3e (threshold %.3e), sustained max|u| %.4e "
          "(threshold %.4e sigma/mu_c)\n",
          res[0].edp, TOL_DP, res[0].usus, TOL_U);
  printf ("  => a 4-cell gap %s the static thresholds at production "
          "property ratios\n", ok4 ? "MEETS" : "DOES NOT MEET");
  printf ("  at the production projection tolerance (%g) the same "
          "configuration gives relative dp error %.3e and sustained max|u| "
          "%.4e: %s\n", TOL_PROD, res[4].edp, res[4].usus,
          (res[4].edp < TOL_DP && res[4].usus < TOL_U) ?
          "also within the static thresholds"
          : "OUTSIDE the static thresholds -- the production projection "
          "tolerance, not the film resolution, is then the binding "
          "constraint");

  /** The gate itself. */

  Result * g = NULL;
  for (int i = 0; i < NRUN; i++)
    if (res[i].gapcells == GAP_GATE && res[i].tol == TOL_BASE)
      g = &res[i];

  bool integrity = true;
  for (int i = 0; i < NRUN; i++)
    if (!res[i].populated || !res[i].pure)
      integrity = false;

  bool pass = (g != NULL) && integrity &&
    (g->edp < TOL_DP) && (g->usus < TOL_U);

  printf ("\nGate (gap = %g cells, TOLERANCE = %g)\n", GAP_GATE, TOL_BASE);
  if (g == NULL)
    printf ("  gated run missing\n");
  else {
    printf ("  relative dp error %.3e < %.3e : %s\n",
            g->edp, TOL_DP, g->edp < TOL_DP ? "yes" : "NO");
    printf ("  sustained max|u| %.4e < %.4e sigma/mu_c : %s\n",
            g->usus, TOL_U, g->usus < TOL_U ? "yes" : "NO");
  }
  printf ("  sampling-region integrity (all runs) : %s\n",
          integrity ? "yes" : "NO");

  /**
  ### Machine-readable output
  */

  FILE * fp = fopen ("staticFilmTube-gap.csv", "w");
  if (fp == NULL) {
    fprintf (stderr, "ERROR: cannot open staticFilmTube-gap.csv\n");
    printf ("FAIL\n");
    exit (1);
  }
  fprintf (fp, "# staticFilmTube.c: static bubble near an embedded tube wall\n");
  fprintf (fp, "# comparator: exact static solution u == 0, dp = 2*sigma/Rd\n");
  fprintf (fp, "# L0=%g\n", Ldom);
  fprintf (fp, "# Rtube=%g\n", (double) Rtube);
  fprintf (fp, "# Xd=%g\n", (double) Xd);
  fprintf (fp, "# N=%d\n", Ngrid);
  fprintf (fp, "# Delta=%.9e\n", delta);
  fprintf (fp, "# sigma=%g\n", f.sigma);
  fprintf (fp, "# rho1=%g\n", rho1);
  fprintf (fp, "# rho2=%g\n", rho2);
  fprintf (fp, "# mu1=%g\n", mu1);
  fprintf (fp, "# mu2=%g\n", mu2);
  fprintf (fp, "# tend=%g\n", (double) tend);
  fprintf (fp, "# DT=%g\n", (double) DTMAX);
  fprintf (fp, "# TOLERANCE_baseline=%g\n", (double) TOL_BASE);
  fprintf (fp, "# TOLERANCE_tight=%g\n", (double) TOL_TIGHT);
  fprintf (fp, "# TOL_DP=%g\n", (double) TOL_DP);
  fprintf (fp, "# TOL_U=%g\n", (double) TOL_U);
  fprintf (fp, "gap_cells,Rd,gap,tolerance,nsteps,dp,dp_exact,rel_error,"
           "umax_tend,umax_tend_cutcells,umax_sustained,umax_peak,t_peak,"
           "umax_mean_2nd_quarter,umax_mean_4th_quarter,umax_ratio,umax_cv,"
           "Rd_eff,dp_exact_eff,n_in,n_out,min_f_in,max_f_out\n");
  for (int i = 0; i < NRUN; i++) {
    Result * r = &res[i];
    fprintf (fp, "%g,%.9e,%.9e,%g,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
             "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%d,%d,%.9e,%.9e\n",
             r->gapcells, r->Rd, r->gap, r->tol, r->nsteps, r->dp, r->dpex,
             r->edp, r->umax, r->umax_cut, r->usus, r->upeak, r->tpeak,
             r->uwin1, r->uwin2, r->uratio, r->ucv, r->Rdeff, r->dpeff,
             r->nin, r->nout, r->fmin_in, r->fmax_out);
  }
  fclose (fp);
  printf ("\nwrote staticFilmTube-gap.csv\n");

  printf ("%s\n", pass ? "PASS" : "FAIL");
  if (!pass)
    exit (1);
  return 0;
}

/**
## Initialisation

The wall is embedded first (which also synchronises the axisymmetric
metric), then the spherical interface is initialised and cleaned inside
the solid. The sphere never reaches the wall, so the cleanup is a no-op
here; it is kept to mirror the production sequence exactly.
*/

event init (t = 0)
{
  tube_solid (Rtube);
  fraction (f, Rd - sqrt (sq(x - Xd) + sq(y)));
  vof_solid_cleanup (f);
}

/**
## Diagnostics

The pressure jump is measured between volume-weighted cell averages of
the two sampling regions defined in the header block. Both regions are
counted and checked for phase purity, so a mis-specified region fails
loudly rather than returning a plausible-looking number. `umax` is taken
over full fluid cells (`cs >= 1`), matching `laplaceEmbedTube.c`;
`umax_cut` additionally includes the cut cells adjacent to the wall,
which is where the height-function stencil is most exposed.

`Rdeff` is the radius recovered from the conserved VOF volume,
$\int f\,\mathrm{d}v = \tfrac{2}{3}R_d^3$ in the axisymmetric measure. It
verifies that the comparator radius is still the radius of the bubble
actually present on the grid.
*/

static void diagnose (Result * r)
{
  double pin = 0., pout = 0., win = 0., wout = 0.;
  double nin = 0., nout = 0., fo = 0., fi = 1.;
  foreach (reduction(+:pin) reduction(+:pout) reduction(+:win)
           reduction(+:wout) reduction(+:nin) reduction(+:nout)
           reduction(max:fo) reduction(min:fi)) {
    if (cs[] >= 1.) {
      double rr = sqrt (sq(x - Xd) + sq(y));
      if (rr < RIN_FRAC*Rd) {
        pin += p[]*dv(); win += dv(); nin += 1.;
        if (f[] < fi) fi = f[];
      }
      else if (x > XOUT) {
        pout += p[]*dv(); wout += dv(); nout += 1.;
        if (f[] > fo) fo = f[];
      }
    }
  }

  double un = 0., unc = 0., vol = 0.;
  foreach (reduction(max:un) reduction(max:unc) reduction(+:vol)) {
    if (cs[] > 0.) {
      double u2 = sqrt (sq(u.x[]) + sq(u.y[]));
      if (u2 > unc) unc = u2;
      if (cs[] >= 1. && u2 > un) un = u2;
      vol += f[]*dv();
    }
  }

  r->Rd = Rd;
  r->gap = Rtube - Rd;
  r->nin = (int) nin;
  r->nout = (int) nout;
  r->fmin_in = fi;
  r->fmax_out = fo;
  r->populated = (nin >= NMIN_IN && nout >= NMIN_OUT);
  r->pure = (nin > 0. && nout > 0. && fi > 1. - FPURE && fo < FPURE);
  r->dp = (win > 0. && wout > 0.) ? pin/win - pout/wout : nan("");
  r->dpex = 2.*f.sigma/Rd;
  r->edp = fabs (r->dp - r->dpex)/r->dpex;
  r->umax = un;
  r->umax_cut = unc;
  r->Rdeff = pow (1.5*vol, 1./3.);
  r->dpeff = 2.*f.sigma/r->Rdeff;
}

/**
The parasitic-current history is stored so that the report can state
whether `max|u|` has plateaued by $t_{end}$ rather than presenting a
snapshot as a converged amplitude. Because the amplitude fluctuates
between samples, the plateau test compares the *mean* over the last
quarter of the run with the mean over the second quarter; the peak over
the second half is reported separately as the conservative amplitude.
*/

event logfile (t += tlog; t <= tend)
{
  Result tmp;
  diagnose (&tmp);
  if (nhist < NHIST) {
    thist[nhist] = t;
    uhist[nhist] = tmp.umax;
    nhist++;
  }
  fprintf (stderr, "%d %g %.6e %.6e %.6e\n",
           run_id, t, tmp.dp, tmp.umax, tmp.umax_cut);
}

event end (t = tend)
{
  Result * r = &res[run_id];
  double gapcells = r->gapcells, tol = r->tol;
  diagnose (r);
  r->gapcells = gapcells;
  r->tol = tol;
  r->nsteps = iter;

  /** Overall peak, the conservative "sustained" amplitude (peak over
  $t \ge t_{end}/2$, the gated statistic), and the windowed means whose
  ratio is the plateau criterion. */

  double upeak = 0., tpeak = 0., usus = 0.;
  double s1 = 0., s2 = 0., sh = 0., shh = 0.;
  int n1 = 0, n2 = 0, nh = 0;
  for (int i = 0; i < nhist; i++) {
    if (uhist[i] > upeak) { upeak = uhist[i]; tpeak = thist[i]; }
    if (thist[i] >= 0.5*tend) {
      if (uhist[i] > usus)
        usus = uhist[i];
      sh += uhist[i]; shh += sq(uhist[i]); nh++;
    }
    if (thist[i] > 0.25*tend && thist[i] <= 0.5*tend) { s1 += uhist[i]; n1++; }
    if (thist[i] > 0.75*tend) { s2 += uhist[i]; n2++; }
  }
  r->upeak = upeak;
  r->tpeak = tpeak;
  r->usus = usus;
  r->uwin1 = (n1 ? s1/n1 : 0.);
  r->uwin2 = (n2 ? s2/n2 : 0.);
  r->uratio = (r->uwin1 > 0. ? r->uwin2/r->uwin1 : HUGE);

  /** Coefficient of variation of `max|u|` over the second half of the
  run. This measures the intrinsic scatter of the statistic and is used
  to calibrate what counts as a *material* change under the
  tolerance-control comparison; without it, any threshold on `max|u|`
  would be arbitrary. */

  if (nh > 1) {
    double mean = sh/nh;
    double var = shh/nh - sq(mean);
    r->ucv = (mean > 0. && var > 0. ? sqrt (var)/mean : 0.);
  }
  else
    r->ucv = 0.;

  fprintf (stderr, "# run %d done: gap %.1f cells, Rd %.6f, TOLERANCE %g, "
           "steps %d, dp %.6e (exact %.6e, rel.err %.3e), max|u| %.4e, "
           "n_in %d, n_out %d\n",
           run_id, r->gapcells, r->Rd, r->tol, r->nsteps, r->dp, r->dpex,
           r->edp, r->umax, r->nin, r->nout);
}
