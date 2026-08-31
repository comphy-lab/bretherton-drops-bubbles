/**
# laplaceEmbedTubeAdapt.c

**Verification** of surface tension inside an embedded tube *on an
adaptive tree*, and a direct test of the three embed/axi/VOF couplings
that `src-local/embed-vof-tube.h` exists to guard.

## Claim

For a static spherical droplet of radius $R_d$ held on the axis of an
embedded tube of radius $R_{tube}$, with matched densities and
viscosities and no gravity, the exact solution of the *implemented*
equations is $\mathbf{u} \equiv 0$ everywhere and a uniform pressure
jump $\Delta p = 2\sigma/R_d$ across the interface. On an adaptive
tree, and provided the grid is re-synchronised after every
`adapt_wavelet()` call, the discrete solution stays within the same
tolerances as the uniform $N = 128$ case.

The comparator is an exact solution of the equations the solver claims
to discretise, so this is *verification*, not validation: it says
nothing about whether those equations describe a real drop in a tube.

## Relation to `laplaceEmbedTube.c`

The physical setup is identical — $L_0 = 2$, $R_d = 0.4$,
$R_{tube} = 0.7$, $X_d = 1$, $\sigma = 1$, $\rho_1 = \rho_2 = 1$,
$\mu_1 = \mu_2 = 1$, `DT` $= 10^{-2}$, `TOLERANCE` $= 10^{-5}$,
$t_{end} = 1$ — and `MAXLEVEL` $= 7$ reproduces that case's finest
spacing $\Delta = L_0/128$. The only difference is the grid: here the
run starts from a coarse `init_grid(1 << MINLEVEL)` and
`adapt_wavelet()` runs every step.

The comparison is therefore direct *at the interface*, which the `f`
criterion pins at `MAXLEVEL`; it is not direct at the wall, which ends
up at `MINLEVEL` — see the refinement criteria below, and the grid
statistics the case prints.

Density and viscosity filtering is deliberately absent: with
$\rho_1 = \rho_2$ and $\mu_1 = \mu_2$ the property functions of
`two-phase-generic.h` are identically unity whether or not the smoothed
field replaces `f`, so `#define FILTERED` would be inert — and the
production solver no longer sets it either.

## Refinement criteria

`adapt_wavelet()` on $\{f, c_s, u_x, u_y\}$ with thresholds
$\{10^{-3}, 10^{-3}, 10^{-2}, 10^{-2}\}$ between `MINLEVEL` $= 4$ and
`MAXLEVEL` $= 7$. The choices:

- `fErr` $= 10^{-3}$ holds the whole droplet interface at `MAXLEVEL`;
- `csErr` $= 10^{-3}$ was intended to hold the tube wall at `MAXLEVEL`.
  **It does not, and this case measures that.** The wall is planar and
  grid-aligned, and `embed.h` sets `cs.prolongation = fraction_refine`,
  which reproduces a planar interface exactly; the wavelet error of
  `cs` is therefore identically zero along a straight tube wall and no
  refinement is requested. Measured here: the `init` event refines the
  wall band to `MAXLEVEL`, the first `adapt_wavelet()` coarsens it
  straight back, and at $t = t_{end}$ **all 16 cut cells sit at
  `MINLEVEL`**. Lowering `csErr` cannot change this — the error is
  zero, not small. This contradicts the remedy suggested in
  `src-local/embed-vof-tube.h` for coupling 3 ("cases must refine the
  wall region (e.g. by adapting on `cs`)"), and it applies equally to
  `simulationCases/bretherton.c`, which uses the same criterion: there,
  wall refinement is delivered by the *interface* criteria once the
  film is close to the wall, not by `cs`. A case that needs the wall
  refined independently of the interface must raise `MINLEVEL` or add
  an explicit geometric criterion.
- `uErr` $= 10^{-2}$ is set an order of magnitude *above* the spurious
  velocity the case is allowed to reach (`TOL_U` $= 10^{-3}$). The
  exact solution has $\mathbf{u} \equiv 0$, so a tighter velocity
  criterion would refine on numerical noise and the mesh would stop
  being a controlled variable. Velocity is kept in the criterion list
  because the production solver adapts on it, but here it is
  deliberately inactive.

The production solver also adapts on the curvature `KAPPA`; that is
omitted here because the interface is a sphere of fixed radius and
`fErr` alone already pins it at `MAXLEVEL`.

Note that `MINLEVEL` $= 4$ means the droplet interior and the far field
coarsen to $\Delta = L_0/16$. Both are regions of uniform pressure, so
the cell-averaged pressures used below are not sensitive to this.

## The three probes

`src-local/embed-vof-tube.h` names three couplings, *all three
triggered by adaptation*. No other case in this repository runs on an
adaptive grid, so until now none of them was exercised in the situation
the header was written for. Each probe below targets one of them and is
evaluated after **every** adaptation step, with the running maximum
carried to the end of the run.

**P1 — metric consistency.** With `AXI` and `EMBED` the metric must
satisfy $c_m = y\,c_s$ and $f_m = y\,f_s$ (with the cut-cell
plane-centre correction of `axi.h`). `axi.h` fills `cm`/`fm` once, in
its `metric` event at `i = 0`, when `cs = 1` everywhere; every later
`solid()` call and every `adapt_wavelet()` — which re-prolongates
`cs`/`fs` — leaves them stale. The probe recomputes reference fields
with the *same* arithmetic as `axi.h`'s `cm_update()` and `fm_update()`
(including the `facet_normal`/`plane_center` correction in cut cells,
the $-\sigma(1 - f_{s,x})$ shift on $x$-faces, and the fact that in
`foreach_face(y)` the coordinate `y` is the face position, not the cell
centre) and reports
$\max|c_m - c_m^{ref}|$ and $\max|f_m - f_m^{ref}|$. The probe iterates
over `foreach()` and `foreach_face()`, i.e. over exactly the domain
`cm_update()`/`fm_update()` write, which includes the faces on the
domain boundary. Both maxima must be at round-off, `TOL_METRIC`
$= 10^{-12}$.

Because the probe recomputes exactly what `embed_axi_metric_sync()`
writes, in the guarded build P1 is close to tautological: it confirms
that the sync was actually invoked and that nothing between the sync
and the probe perturbed the metric. Its discriminating power is
demonstrated by the negative control, where the same probe fires.

**P2 — no dispersed phase inside the solid.** `vof.h` only updates
cells with `cs > 0`, but adaptation prolongates `f` with no knowledge
of `cs`, so refined solid cells can inherit non-zero `f` from fluid
neighbours, corrupting the (non-embed-aware) height-function columns
near the wall. The probe reports $\max|f|$ over cells with
$c_s \le 0$ and requires it to be exactly zero.

**P3 — Young–Laplace.** Relative error of the measured pressure jump
below `TOL_DP` $= 2\times 10^{-2}$ and $\max|\mathbf{u}|$ below
`TOL_U` $= 10^{-3}$ at $t = t_{end}$, using the same thresholds and the
same measurement regions as the uniform case: cell averages of `p` over
$r < R_d/2$ inside and over $r > 3R_d/2$, $y < R_{tube} - 0.1$, full
fluid cells only, outside. This is the physics-level consequence of
couplings 1–3 being handled correctly.

## Negative control

Compile with `-DSKIP_EMBED_GUARDS` to omit the post-adaptation
`embed_axi_metric_sync()` and `vof_solid_cleanup()` calls — and only
those; the `init` event still synchronises through `tube_solid()`, so
the initial state is identical in both builds and the control isolates
the *adaptation* coupling.

~~~bash
qcc -I../src-local -O2 -Wall -disable-dimensions \
    laplaceEmbedTubeAdapt.c -o laplaceEmbedTubeAdapt -lm
qcc -I../src-local -O2 -Wall -disable-dimensions -DSKIP_EMBED_GUARDS \
    laplaceEmbedTubeAdapt.c -o laplaceEmbedTubeAdaptNoGuard -lm
~~~

The guard-less build is expected to fail, and to fail *on P1 or P2*,
i.e. for the reason the header states, rather than merely by producing
NaN. The probes therefore run before any solver step consumes the
corrupted state, and the first firing iteration and magnitude are
recorded and printed. If the guard-less build were to pass, that would
be a finding about the header — not a reason to weaken this case.

Observed behaviour (this machine, serial, Basilisk as of writing):

- the guard-less build **fails, on P1**, at the very first adaptation
  (`i = 0`), with $\max|f_m - f_m^{ref}| = 6.25\times10^{-4}$ against a
  local $f_m \approx 0.28$, i.e. a relative metric error of about
  $2\times10^{-3}$. $\max|c_m - c_m^{ref}|$ stays exactly zero: `cm`
  has a dedicated axisymmetric prolongation (`refine_cm_axi`) and its
  restriction happens to be exact for a planar wall, whereas the
  $x$-face metric is not restored;
- the violation persists, unchanged in magnitude, for all 921 steps;
- **P2 never fires** (see the limits below);
- **P3 is bit-identical between the two builds**: same $\Delta p$, same
  $\max|\mathbf{u}|$ to all printed digits, same step count, same leaf
  count. Instrumenting the offending face shows why: in this geometry
  the only stale face is the one on the outflow boundary $x = L_0$ in
  the wall cut-cell row, where the default no-flux condition makes the
  face flux zero regardless of $f_m$.

So the negative control does fail for the reason the header claims —
the metric is left inconsistent by adaptation — but in *this*
configuration that inconsistency has no measurable dynamical
consequence. The case therefore establishes that
`embed_axi_metric_sync()` is doing something real and necessary for
metric consistency; it does **not** establish that omitting it would
corrupt a production Bretherton run. That would need a configuration in
which stale interior faces carry non-zero flux.

## Outputs

- the running maxima of P1 and P2 and their first-firing iterations;
- $\Delta p$, its exact value, the relative error and
  $\max|\mathbf{u}|$ at $t = t_{end}$;
- final grid statistics: leaf count, level range, finest $\Delta$, and
  the number and level range of the cut cells — the last of these is
  what exposes the `cs`-refinement behaviour discussed above;
- a time series of $\Delta p$, $\max|\mathbf{u}|$ and the P1/P2 maxima
  on stderr;
- a bare `PASS` or `FAIL` line on stdout, with `exit(1)` on failure.

## Limits

- One geometry, one droplet radius, one tube radius, one
  interface–wall separation $R_{tube} - R_d = 0.3$: about 19 cells at
  `MAXLEVEL` on the interface side, about 2.4 cells at `MINLEVEL` on
  the wall side. Coupling 3 (interface–wall separation) is therefore
  *satisfied by construction* here, not tested: this case establishes
  nothing about what happens when the interface approaches the wall to
  within a few cells, which is exactly the Bretherton film regime.
- The wall ends up at `MINLEVEL` (see the refinement criteria). That
  is a measured property of `adapt_wavelet` on `cs` for a planar
  boundary, reported rather than worked around, but it does mean the
  cut cells here are coarse and that the metric coupling is exercised
  only through the single 7 → 4 coarsening of the wall band during the
  first adaptation steps, not through sustained churn of the cut cells.
- For the same reason the dispersed phase never reaches the wall, so
  `f` is zero in the solid region throughout and **P2 has no
  opportunity to fire**. It is a necessary condition that the guarded
  build satisfies, and a genuine regression guard, but the negative
  control here can only exercise P1. A case in which the drop touches
  or nearly touches the wall would be needed to exercise P2.
- Matched properties only ($\rho_1 = \rho_2$, $\mu_1 = \mu_2$):
  nothing is established about the density and viscosity contrasts the
  production solver runs at.
- One resolution: `MAXLEVEL` $= 7$ only, so this is a tolerance check
  at fixed resolution, not a convergence study. The uniform
  refinement sequence lives in `laplaceEmbedTube.c`.
- The measured jump is a difference of cell averages over finite
  regions, not a pointwise jump, so it carries its own discretisation
  error.
- P1 is a consistency check of the *live* metric against `axi.h`'s own
  formulae. It does not verify that those formulae are themselves the
  correct axisymmetric metric; that is `axi.h`'s claim, guarded
  upstream by
  [`src/test/missing_metric.c`](http://basilisk.fr/src/test/missing_metric.c).
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
#define tend 1.0

#define MAXLEVEL 7
#define MINLEVEL 4

#define fErr 1e-3
#define csErr 1e-3
#define uErr 1e-2

#define TOL_DP 2e-2
#define TOL_U 1e-3
#define TOL_METRIC 1e-12

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

/**
## Probe state

Running maxima and the iteration at which each probe first exceeds its
threshold (`-1` while it has never fired).
*/
static double dcmMax = 0., dfmMax = 0., fSolidMax = 0.;
static int iFirstMetric = -1, iFirstSolid = -1;
static double vFirstMetric = 0., vFirstSolid = 0.;

int main()
{
  L0 = 2.0;
  init_grid (1 << MINLEVEL);

  rho1 = 1.; mu1 = 1.;
  rho2 = 1.; mu2 = 1.;
  f.sigma = 1.;

  TOLERANCE = 1e-5;
  DT = 1e-2;

  run();
}

/**
## Initialisation

Refine the interface and the wall to `MAXLEVEL` *before* embedding the
solid and initialising `f`, so that both are laid down on the mesh they
will be carried on. `tube_solid()` synchronises the metric; this
happens in both builds, so the two differ only in what they do after
`adapt_wavelet()`.
*/
event init (t = 0)
{
  refine (fabs (sqrt (sq(x - Xd) + sq(y)) - Rd) < 0.1 && level < MAXLEVEL);
  refine (fabs (y - Rtube) < 0.1 && level < MAXLEVEL);
  tube_solid (Rtube);
  fraction (f, Rd - sqrt (sq(x - Xd) + sq(y)));
  vof_solid_cleanup (f);
}

/**
## P1 and P2: the coupling probes

`metric_deviation()` recomputes `cm` and `fm` from the live `cs`/`fs`
with the same arithmetic as `cm_update()`/`fm_update()` in
`$BASILISK/axi.h` and returns the maximum absolute deviation of the
live fields from that reference. `cm` and `fm` are read through
plain-field aliases, following `axi.h`'s own idiom for these
`(const)`-declared metric fields.

`solid_fraction_max()` returns $\max|f|$ over full-solid cells: the
absolute value, so that a *negative* fraction prolongated into the
solid is caught as well as a positive one.
*/
static void metric_deviation (double * dcm, double * dfm)
{
  double ecm = 0., efm = 0.;
#if defined(AXI) && defined(EMBED)
  scalar cmref[];
  face vector fmref[];

  foreach() {
    if (cs[] > 0. && cs[] < 1.) {
      coord p, n = facet_normal (point, cs, fs);
      double alpha = plane_alpha (cs[], n);
      plane_center (n, alpha, cs[], &p);
      cmref[] = (y + Delta*p.y)*cs[];
    }
    else
      cmref[] = y*cs[];
  }

  /**
  In `foreach_face(x)` the coordinate `y` is the cell-centre radius (an
  $x$-face has the radius of its cell); in `foreach_face(y)` it is the
  face radius, which vanishes on the axis — hence the `max(y, 1e-20)`
  floor, copied from `fm_update()`. */

  foreach_face(x) {
    double sig = 0.;
    if (cs[] > 0. && cs[] < 1.) {
      coord n = facet_normal (point, cs, fs);
      sig = sign(n.y)*Delta/2.;
    }
    fmref.x[] = (y - sig*(1. - fs.x[]))*fs.x[];
  }
  foreach_face(y)
    fmref.y[] = max(y, 1e-20)*fs.y[];

  scalar cmv = cm;
  face vector fmv = fm;

  foreach (reduction(max:ecm))
    if (fabs (cmv[] - cmref[]) > ecm)
      ecm = fabs (cmv[] - cmref[]);
  foreach_face (x, reduction(max:efm))
    if (fabs (fmv.x[] - fmref.x[]) > efm)
      efm = fabs (fmv.x[] - fmref.x[]);
  foreach_face (y, reduction(max:efm))
    if (fabs (fmv.y[] - fmref.y[]) > efm)
      efm = fabs (fmv.y[] - fmref.y[]);
#endif
  *dcm = ecm; *dfm = efm;
}

static double solid_fraction_max (void)
{
  double fm2 = 0.;
  foreach (reduction(max:fm2))
    if (cs[] <= 0. && fabs (f[]) > fm2)
      fm2 = fabs (f[]);
  return fm2;
}

/**
Updates the running maxima and records the first violation of each
predicate. Called after every adaptation.
*/
static void probe_couplings (int iter)
{
  double dcm, dfm;
  metric_deviation (&dcm, &dfm);
  double fsol = solid_fraction_max();

  if (dcm > dcmMax) dcmMax = dcm;
  if (dfm > dfmMax) dfmMax = dfm;
  if (fsol > fSolidMax) fSolidMax = fsol;

  if (iFirstMetric < 0 && (dcm >= TOL_METRIC || dfm >= TOL_METRIC)) {
    iFirstMetric = iter;
    vFirstMetric = max (dcm, dfm);
    fprintf (ferr, "# P1 fired at i = %d: max|cm - y*cs| = %.6e, "
             "max|fm - y*fs| = %.6e\n", iter, dcm, dfm);
  }
  if (iFirstSolid < 0 && fsol > 0.) {
    iFirstSolid = iter;
    vFirstSolid = fsol;
    fprintf (ferr, "# P2 fired at i = %d: max|f| over cs <= 0 = %.6e\n",
             iter, fsol);
  }
}

/**
## P3: pressure jump and spurious currents

Identical to `laplaceEmbedTube.c`: cell averages of `p` well inside
($r < R_d/2$) and well outside ($r > 3R_d/2$, away from the wall) the
droplet, and the maximum velocity magnitude over every cell with
$c_s > 0$, cut cells included.
*/
static void dp_umax (double * dp, double * umax)
{
  double pin = 0., pout = 0., win = 0., wout = 0., un = 0.;
  foreach (reduction(+:pin) reduction(+:pout)
           reduction(+:win) reduction(+:wout) reduction(max:un)) {

    /**
    The spurious-current maximum must include cut cells. Those are the
    cells adjacent to the embedded wall, and they are where the
    height-function curvature error this case measures is largest;
    restricting `un` to full cells would let a near-wall parasitic
    current exceed `TOL_U` unseen. The pressure averages stay on pure
    fluid cells, where a cell-average of `p` is meaningful. */

    if (cs[] > 0.) {
      double u2 = sqrt (sq(u.x[]) + sq(u.y[]));
      if (u2 > un)
        un = u2;
    }
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
    }
  }
  *dp = (win > 0. && wout > 0.) ? pin/win - pout/wout : nodata;
  *umax = un;
}

/**
## Adaptation

The whole point of the case. `adapt_wavelet()` re-prolongates `cs`,
`fs` and `f` with no knowledge of the axisymmetric metric or of which
cells are solid; `embed_axi_metric_sync()` and `vof_solid_cleanup()`
repair both, exactly as `simulationCases/bretherton.c` does. The
`SKIP_EMBED_GUARDS` build omits them and nothing else.
*/
event adapt (i++)
{
  adapt_wavelet ((scalar *){f, cs, u.x, u.y},
                 (double[]){fErr, csErr, uErr, uErr},
                 MAXLEVEL, MINLEVEL);
#if !defined(SKIP_EMBED_GUARDS)
  embed_axi_metric_sync();
  vof_solid_cleanup (f);
#endif
  probe_couplings (i);
}

/**
## Reporting

`summarise()` prints every measured quantity, the verdict for each
probe and a bare `PASS`/`FAIL` line, then exits with status 1 on
failure. It is called at $t = t_{end}$, and early if the solution has
become non-finite — in that case the recorded first-firing data are
still printed, which is what distinguishes "failed on P1/P2" from
"merely blew up".
*/
static void grid_stats (long * ncells, int * lmin, int * lmax)
{
  long n = 0;
  int lo = 30, hi = 0;
  foreach (reduction(+:n) reduction(min:lo) reduction(max:hi)) {
    n++;
    if (level < lo) lo = level;
    if (level > hi) hi = level;
  }
  *ncells = n; *lmin = lo; *lmax = hi;
}

/**
The cut cells carry the embedded wall, so their count and level range
say directly whether the wall is resolved — which is what coupling 3
of `src-local/embed-vof-tube.h` asks a case to monitor.
*/
static void cut_cell_stats (long * ncut, int * lmin, int * lmax)
{
  long n = 0;
  int lo = 30, hi = 0;
  foreach (reduction(+:n) reduction(min:lo) reduction(max:hi))
    if (cs[] > 0. && cs[] < 1.) {
      n++;
      if (level < lo) lo = level;
      if (level > hi) hi = level;
    }
  *ncut = n; *lmin = lo; *lmax = hi;
}

static void summarise (const char * note, int iter, bool fatal)
{
  double dp, umax;
  dp_umax (&dp, &umax);
  double dpex = 2./Rd;
  double edp = fabs (dp - dpex)/dpex;

  long ncells, ncut; int lmin, lmax, lcmin, lcmax;
  grid_stats (&ncells, &lmin, &lmax);
  cut_cell_stats (&ncut, &lcmin, &lcmax);

  bool p1 = (dcmMax < TOL_METRIC && dfmMax < TOL_METRIC);
  bool p2 = (fSolidMax == 0.);
  bool p3 = (edp < TOL_DP && umax < TOL_U);

#if defined(SKIP_EMBED_GUARDS)
  printf ("build: SKIP_EMBED_GUARDS (negative control, "
          "no post-adaptation sync/cleanup)\n");
#else
  printf ("build: guarded (embed_axi_metric_sync + vof_solid_cleanup "
          "after every adapt_wavelet)\n");
#endif
  printf ("stopped at: %s (t = %g, i = %d)\n", note, t, iter);
  printf ("grid: %ld leaves, levels %d-%d, finest Delta = %.6e\n",
          ncells, lmin, lmax, L0/(1 << lmax));
  if (ncut > 0)
    printf ("cut cells (embedded wall): %ld, levels %d-%d, "
            "Delta at wall = %.6e\n", ncut, lcmin, lcmax,
            L0/(1 << lcmax));
  else
    printf ("cut cells (embedded wall): none\n");

  printf ("P1 metric consistency: max|cm - cm_ref| = %.6e, "
          "max|fm - fm_ref| = %.6e (tol %.1e) -> %s\n",
          dcmMax, dfmMax, (double) TOL_METRIC, p1 ? "pass" : "FAIL");
  if (iFirstMetric >= 0)
    printf ("  first violated at i = %d, magnitude %.6e\n",
            iFirstMetric, vFirstMetric);
  printf ("P2 dispersed phase in solid: max|f| over cs <= 0 = %.6e "
          "(required exactly 0) -> %s\n", fSolidMax, p2 ? "pass" : "FAIL");
  if (iFirstSolid >= 0)
    printf ("  first violated at i = %d, magnitude %.6e\n",
            iFirstSolid, vFirstSolid);
  printf ("P3 Young-Laplace: dp = %.6e (exact %.6e, rel. error %.3e, "
          "tol %.1e)\n", dp, dpex, edp, (double) TOL_DP);
  printf ("   max spurious velocity = %.6e at TOLERANCE %g (tol %.1e) -> %s\n",
          umax, TOLERANCE, (double) TOL_U, p3 ? "pass" : "FAIL");

  /**
  A non-finite state is a failure in its own right, even with every probe
  green: the probes sample fluid cells, so NaNs confined to solid or cut
  cells would otherwise be reported as a pass. */

  bool ok = p1 && p2 && p3 && !fatal;
  if (fatal)
    printf ("divergence: the solution is not finite, so the probes above "
            "are not meaningful\n");
  printf (ok ? "PASS\n" : "FAIL\n");
  fflush (stdout);
  if (!ok)
    exit (1);
}

/**
Divergence guard: if the solution is no longer finite there is nothing
left to measure, so report what was recorded and stop. This never
triggers in the guarded build.
*/
static bool nonfinite_state (void)
{
  double bad = 0.;
  foreach (reduction(max:bad))
    if (!isfinite (u.x[]) || !isfinite (u.y[]) || !isfinite (p[]))
      bad = 1.;
  return bad > 0.;
}

event logfile (t += 0.1; t <= tend)
{
  double dp, umax;
  dp_umax (&dp, &umax);
  fprintf (stderr, "%g %g %.6e %.6e %.6e %.6e %.6e\n",
           TOLERANCE, t, dp, umax, dcmMax, dfmMax, fSolidMax);
}

event sanity (i++)
{
  if (nonfinite_state())
    summarise ("non-finite solution", i, true);
}

event end (t = tend)
{
  summarise ("t = tend", i, false);
}
