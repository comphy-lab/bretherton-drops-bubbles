/**
# laplaceEmbedTube.c

**Verification** of surface tension inside an embedded tube, as a
grid-refinement sequence.

## Claim

For a static spherical droplet of radius $R_d$ held on the axis of an
embedded tube of radius $R_{tube}$, with matched densities and
viscosities and no gravity, the exact solution of the *implemented*
equations is $\mathbf{u} \equiv 0$ everywhere and a uniform pressure
jump $\Delta p = 2\sigma/R_d$ across the interface. This case measures
how the discrete solution approaches that comparator under uniform
refinement.

The comparator is an exact solution of the equations the solver claims
to discretise, so this is *verification*, not validation: it says
nothing about whether those equations describe a real drop in a tube.

Any velocity at $t = t_{end}$ is a spurious (parasitic) current; any
deviation of the measured jump from $2\sigma/R_d$ measures the
curvature error of the height-function method operating near an
embedded boundary it knows nothing about.

This exercises the production stack of `simulationCases/bretherton.c`
— embed + axi + centered + two-phase + tension — in a configuration
with a known solution. The droplet surface stays $R_{tube} - R_d$ away
from the wall, mimicking the wetting-film separation of the Bretherton
configuration.

## Refinement sequence

Three uniform levels, $N = 64, 128, 256$, driven the idiomatic Basilisk
way by `for (N = 64; N <= 256; N *= 2) run();` in `main()` (see e.g.
`$BASILISK/test/couette.c`, `$BASILISK/test/spurious.c`). `run()` calls
`init_grid(N)` on entry and `free_grid()` on exit, so every field is
reallocated and zeroed between levels; only the `main()`-set parameters
and the file-scope embedded-boundary conditions — installed once in the
generated `_init_solver()`, not per run — survive. This was checked
against a stand-alone $N = 128$ run: the swept and stand-alone
$N = 128$ results agree to all printed digits.

The recorded sequence gives relative pressure-jump errors of
$2.74\times 10^{-3}$, $6.05\times 10^{-4}$ and $1.92\times 10^{-4}$,
i.e. observed orders $p = 2.18$ then $1.66$. The rate is not constant
over these two intervals and no mechanism for that is established
here; three levels are too few to identify the asymptotic order.

The physical setup is identical at every level: $L_0 = 2$,
$R_d = 0.4$, $R_{tube} = 0.7$, $\sigma = 1$, $\rho_1 = \rho_2 = 1$,
$\mu_1 = \mu_2 = 1$, `DT` $= 10^{-2}$, `TOLERANCE` $= 10^{-5}$,
$t_{end} = 1$. Only $\Delta = L_0/N$ changes.

Density and viscosity filtering is deliberately absent here: with
$\rho_1 = \rho_2$ and $\mu_1 = \mu_2$, `two-phase-generic.h` evaluates
`rho(f) = clamp(f)*(rho1 - rho2) + rho2` and `mu(f)` likewise, both
identically unity whether or not the smoothed field `sf` replaces `f`,
so `#define FILTERED` would be inert — and the production solver no
longer sets it either.

## Outputs

- one line per level on stdout: $N$, $\Delta$, cells across $R_d$,
  measured $\Delta p$, exact $\Delta p$, relative error, and
  $\max|\mathbf{u}|$ at $t = t_{end}$;
- the observed order of convergence of the pressure-jump error between
  successive levels, $p = \log_2(e_{coarse}/e_{fine})$;
- `laplaceEmbedTube-refinement.csv` in the working directory: a
  commented metadata header, a column header row, then one machine-
  readable row per level;
- the time series of $\Delta p$ and $\max|\mathbf{u}|$ on stderr,
  prefixed by $N$.

## Pass predicate

1. at the finest level, relative $\Delta p$ error $< 2\%$ **and**
   $\max|\mathbf{u}| < 10^{-3}\,\sigma/\mu$; **and**
2. the relative $\Delta p$ error decreases monotonically with
   refinement.

The spurious velocity is deliberately **not** gated on decreasing.
Parasitic currents in a static drop frequently plateau or grow under
refinement, because the height-function curvature error and the
discrete force/pressure balance do not refine at the same rate, so
making the case fail on a growing $\max|\mathbf{u}|$ would be a
statement about the test rather than about the solver. It is measured,
recorded and printed at every level, but its trend does not decide
pass or fail.

At the three levels run here it does in fact decrease monotonically —
$1.82\times 10^{-4}$, $4.14\times 10^{-5}$, $1.29\times 10^{-5}$ at
$N = 64, 128, 256$, i.e. rates $p \simeq 2.13$ then $1.68$ — measured
at the fixed finite time $t = 1$, not at a steady state. That is
consistent with converging parasitic currents but does not establish
it: two rate estimates over a factor of four in $\Delta$, at a time
where the currents are still decaying (see the stderr time series),
cannot distinguish decay of the currents from decay of the transient.
The case therefore bounds $\max|\mathbf{u}|$ below $10^{-3}$ at the
resolutions tested and reports its trend; it does not certify
convergence of the spurious currents. Should that trend reverse at
higher resolution, the case will report it and still pass.

## Limits

- One geometry, one droplet radius, one tube radius, one interface–wall
  separation. Nothing is established about other $R_d/R_{tube}$.
- Uniform grid only: the adaptive path (`adapt_wavelet()` plus the
  `embed_axi_metric_sync()` / `vof_solid_cleanup()` calls it requires)
  is *not* exercised here — see `laplaceEmbedTubeAdapt.c`.
- Matched properties only: nothing is established about the
  density- and viscosity-contrast cases the production solver runs.
- Three levels give two order estimates; that is enough to see a trend,
  not enough to claim an asymptotic rate.
- The measured jump is a difference of cell averages over finite
  regions, not a pointwise jump, so it carries its own
  discretisation error.
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "embed-vof-tube.h"

#define Rtube 0.7
#define Rd 0.4
#define Xd 1.0
#define SIGMA 1.
#define tend 1.0
#define TOL_DP 2e-2
#define TOL_U 1e-3

/**
The refinement sequence: `NLEVELS` uniform grids from `NCOARSE` to
`NFINE`, doubling. */

#define NCOARSE 64
#define NFINE 256
#define NLEVELS 3

int lev_N[NLEVELS];
double lev_dp[NLEVELS], lev_umax[NLEVELS];
int ilev = 0;

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

int main()
{
  L0 = 2.0;

  rho1 = 1.; mu1 = 1.;
  rho2 = 1.; mu2 = 1.;
  f.sigma = SIGMA;

  TOLERANCE = 1e-5;
  DT = 1e-2;

  /**
  Each `run()` reallocates the grid, so the fields are reinitialised
  between levels; the `init` event below rebuilds the solid and the
  interface from scratch every time. */

  for (N = NCOARSE; N <= NFINE; N *= 2)
    run();

  /**
  ## Post-processing of the sequence

  Everything below is plain C on the recorded per-level scalars; the
  grid is already freed. */

  double dpex = 2.*SIGMA/Rd;
  double e[NLEVELS];

  if (ilev != NLEVELS) {
    fprintf (stderr, "ERROR: recorded %d levels, expected %d\n",
             ilev, NLEVELS);
    printf ("FAIL\n");
    exit (1);
  }

  for (int i = 0; i < NLEVELS; i++)
    e[i] = fabs (lev_dp[i] - dpex)/dpex;

  /**
  ### Machine-readable refinement table */

  FILE * fp = fopen ("laplaceEmbedTube-refinement.csv", "w");
  if (!fp) {
    fprintf (stderr, "ERROR: cannot open laplaceEmbedTube-refinement.csv\n");
    printf ("FAIL\n");
    exit (1);
  }
  fprintf (fp, "# case: laplaceEmbedTube\n");
  fprintf (fp, "# static spherical droplet on the axis of an embedded tube\n");
  fprintf (fp, "# exact solution: u = 0, dp = 2*sigma/Rd = %.17g\n", dpex);
  fprintf (fp, "# L0=%.17g Rd=%.17g Rtube=%.17g sigma=%.17g\n",
           L0, (double) Rd, (double) Rtube, (double) SIGMA);
  fprintf (fp, "# rho1=%.17g rho2=%.17g mu1=%.17g mu2=%.17g\n",
           rho1, rho2, mu1, mu2);
  fprintf (fp, "# tend=%.17g DT=%.17g TOLERANCE=%.17g\n",
           (double) tend, DT, TOLERANCE);
  fprintf (fp, "# TOL_DP=%.17g TOL_U=%.17g\n",
           (double) TOL_DP, (double) TOL_U);
  fprintf (fp, "N,Delta,cells_across_Rd,dp,dp_exact,rel_error,max_u\n");
  for (int i = 0; i < NLEVELS; i++) {
    double Delta_i = L0/lev_N[i];
    fprintf (fp, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
             lev_N[i], Delta_i, Rd/Delta_i, lev_dp[i], dpex, e[i],
             lev_umax[i]);
  }
  fclose (fp);

  /**
  ### Human-readable summary */

  printf ("\n## laplaceEmbedTube refinement sequence"
          " (exact dp = %.6e)\n", dpex);
  printf ("%6s %10s %10s %14s %12s %12s\n",
          "N", "Delta", "Rd/Delta", "dp", "rel. error", "max|u|");
  for (int i = 0; i < NLEVELS; i++) {
    double Delta_i = L0/lev_N[i];
    printf ("%6d %10.6f %10.3f %14.6e %12.3e %12.3e\n",
            lev_N[i], Delta_i, Rd/Delta_i, lev_dp[i], e[i], lev_umax[i]);
  }

  /**
  ### Observed order of convergence of the pressure-jump error */

  for (int i = 1; i < NLEVELS; i++) {
    if (e[i] > 0. && e[i-1] > 0.)
      printf ("observed order, N = %d -> %d: p = %.3f\n",
              lev_N[i-1], lev_N[i], log2 (e[i-1]/e[i]));
    else
      printf ("observed order, N = %d -> %d: undefined"
              " (vanishing error)\n", lev_N[i-1], lev_N[i]);
  }

  /**
  ### Predicates

  The pressure-jump error must decrease monotonically; the spurious
  velocity is reported but not gated, because parasitic currents need
  not decrease under refinement. */

  bool monotone = true;
  for (int i = 1; i < NLEVELS; i++)
    if (!(e[i] < e[i-1]))
      monotone = false;

  bool umax_decreasing = true;
  for (int i = 1; i < NLEVELS; i++)
    if (!(lev_umax[i] < lev_umax[i-1]))
      umax_decreasing = false;

  int last = NLEVELS - 1;
  bool ok_dp = e[last] < TOL_DP;
  bool ok_u = lev_umax[last] < TOL_U;
  bool pass = ok_dp && ok_u && monotone;

  printf ("finest level N = %d: rel. dp error %.3e (tol %.1e) -> %s\n",
          lev_N[last], e[last], (double) TOL_DP, ok_dp ? "ok" : "EXCEEDED");
  printf ("finest level N = %d: max|u| %.3e (tol %.1e) -> %s\n",
          lev_N[last], lev_umax[last], (double) TOL_U,
          ok_u ? "ok" : "EXCEEDED");
  printf ("dp error monotonically decreasing: %s\n",
          monotone ? "yes" : "NO");
  printf ("max|u| monotonically decreasing: %s"
          " (reported only, not part of the pass predicate)\n",
          umax_decreasing ? "yes" : "no");
  if (!umax_decreasing)
    printf ("note: this case does not demonstrate convergence of the"
            " spurious currents.\n");

  printf (pass ? "PASS\n" : "FAIL\n");
  if (!pass)
    exit (1);

  return 0;
}

event init (t = 0)
{
  tube_solid (Rtube);
  fraction (f, Rd - sqrt (sq(x - Xd) + sq(y)));
  vof_solid_cleanup (f);
}

/**
## Diagnostics

Spurious-current amplitude and pressure jump, logged over time and
recorded at $t = t_{end}$. The pressure jump is measured between cell
averages well inside ($r < R_d/2$) and well outside ($r > 3R_d/2$,
fluid only) the droplet.
*/
static void dp_umax (double * dp, double * umax)
{
  double pin = 0., pout = 0., win = 0., wout = 0., un = 0.;
  foreach (reduction(+:pin) reduction(+:pout)
           reduction(+:win) reduction(+:wout) reduction(max:un)) {
    if (cs[] >= 1.) {
      double r = sqrt (sq(x - Xd) + sq(y));
      if (r < 0.5*Rd) {
        pin += p[]*dv();
        win += dv();
      }
      else if (r > 1.5*Rd && y < Rtube - 0.1) {
        pout += p[]*dv();
        wout += dv();
      }
      double u2 = sqrt (sq(u.x[]) + sq(u.y[]));
      if (u2 > un)
        un = u2;
    }
  }
  *dp = pin/win - pout/wout;
  *umax = un;
}

event logfile (t += 0.1; t <= tend)
{
  double dp, umax;
  dp_umax (&dp, &umax);
  fprintf (stderr, "%d %g %.6e %.6e\n", N, t, dp, umax);
}

/**
At the end of each level we record the raw measurements; no pass/fail
decision is taken here, so the sequence always runs to completion and
the CSV always contains every level.
*/
event end (t = tend)
{
  double dp, umax;
  dp_umax (&dp, &umax);
  if (ilev < NLEVELS) {
    lev_N[ilev] = N;
    lev_dp[ilev] = dp;
    lev_umax[ilev] = umax;
    ilev++;
  }
  printf ("N = %d: pressure jump %.6e (exact %.6e), max|u| %.6e\n",
          N, dp, 2.*SIGMA/Rd, umax);
}
