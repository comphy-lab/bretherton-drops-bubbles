/**
# bretherton.c

Axisymmetric two-phase simulation of a long drop or bubble translating in
a liquid-filled capillary tube, after
[Bretherton (1961)](https://doi.org/10.1017/S0022112061000160).

The tube wall is an *embedded boundary* (`embed.h`), fully wetted by the
continuous phase: the drop/bubble never touches the wall and rides on a
thin lubricating film. The interface is tracked with VOF (`two-phase.h`).
The embed/VOF couplings that Basilisk does not handle automatically are
collected in `src-local/embed-vof-tube.h`.

## Non-dimensionalisation

Repeating variables: surface tension $\sigma$, continuous-phase dynamic
viscosity $\mu_c$, and the volume-equivalent drop/bubble radius
$R = (3V/4\pi)^{1/3}$. Hence

- length scale $R$,
- velocity scale $V_\mu = \sigma/\mu_c$ (visco-capillary velocity),
- time scale $\tau = \mu_c R/\sigma$,
- pressure scale $\sigma/R$.

In these units the *dimensionless bubble tip velocity is itself the
capillary number* $Ca_b = \mu_c U_b/\sigma$.

Control parameters (all runtime keys):

- `Ca`: imposed mean inlet velocity $\mu_c U/\sigma$ (capillary number of
  the driving flow),
- `La`: Laplace number $\rho_c \sigma R/\mu_c^2$ (sets the continuous
  phase density; Bretherton's analysis assumes the visco-capillary limit,
  so keep $La\,Ca \ll 1$),
- `muR`: viscosity ratio $\mu_d/\mu_c$ (bubble: $10^{-2}$; drop: $\geq 1$),
- `rhoR`: density ratio $\rho_d/\rho_c$ (bubble: $10^{-3}$; drop: 1),
- `Rtube`: tube radius in units of $R$ (must be $< 1$ for a confined,
  elongated drop/bubble),
- `Rb0frac`: initial capsule radius as a fraction of `Rtube`,
- `xRear`: initial distance of the rear meniscus tip from the inlet.

Validation targets (small $Ca_b$): film thickness
$b/R_{tube} \simeq 1.34\,Ca_b^{2/3}$ and speed excess
$W = (U_b - U)/U_b \simeq 1.29\,(3 Ca_b)^{2/3}$
[Bretherton (1961), eq. 2]; at moderate $Ca_b$ compare with the
Aussillous & Quéré (2000) fit
$b/R_{tube} = 1.34\,Ca_b^{2/3}/(1 + 2.5\cdot1.34\,Ca_b^{2/3})$.

## Input parameters

Runtime keys via `src-local/params.h` (defaults in brackets): `CaseNo`
[1000], `MAXlevel` [10], `MINlevel` [4], `Ca` [0.05], `La` [1], `muR`
[0.01], `rhoR` [0.001], `Rtube` [0.7], `Rb0frac` [0.8], `xRear` [1.0],
`Ldomain` [16], `travelR` [10], `tmax` [`travelR`/`Ca`],
`tsnap` [`tmax`/200], `tRamp` [1], `dtmax` [0.01], `bTol` [2e-3],
`advWin` [0.25], `advMin` [1.0], `convHold` [3], `uRel` [1e-2],
`dRel` [1e-2], `csErr` [1e-2].

`tmax` and `tsnap` derive from the capillary number by default. The
bubble advances at $U = Ca$, so a fixed run time gives a different
travel distance at every $Ca$, and the film only reaches its steady
thickness after several bubble lengths of advance. `travelR` sets that
advance in units of $R$; the domain holds about 10.4. In practice the
run normally ends earlier, when the film thickness stops changing: see
the convergence stop in `logWriting`. `tmax` is then a cap rather than
the expected duration.
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "params.h"
#include "embed-vof-tube.h"

/**
## Adaptivity controls
*/
#define fErr (1e-3)   // error tolerance in f VOF
#define KErr (1e-4)   // error tolerance in VOF curvature

/**
The velocity and strain-rate tolerances are *relative to the imposed
scales*, not absolute. The velocity scale here is $U = Ca$, which spans
0.002 to 0.05 across the campaign, so a fixed tolerance inherited from
problems where $U \sim 1$ is inert at the low-$Ca$ end: at $Ca = 0.005$
an absolute `1e-2` exceeds 45% of the peak speed, and the mesh follows
the interface alone while the bulk velocity and dissipation fields stay
coarse. That under-resolution grows monotonically as $Ca$ falls.

`uRel` is the velocity tolerance as a fraction of $U$, and `dRel` the
strain-rate tolerance as a fraction of the bulk shear rate
$U/R_{tube}$. */

/**
## Global runtime variables
*/
int MAXlevel, MINlevel, CaseNo;
double Ca, La, muR, rhoR, Rtube, Rb0frac, xRear, Ldomain;

/**
`tsnap` must be non-zero *statically*: Basilisk classifies event
expressions (`t += tsnap` vs conditions) before `main()` assigns the
runtime parameters, and a zero increment is misread as a second
condition. */

double tmax = 200., tsnap = 1., tRamp = 1., travelR = 10.;
double Rb0, Lcyl, Xb0, vol0;
double bTol, advWin, advMin, uRel, dRel, VelErr, DErr, csErr;
int convHold;

char nameOut[128], dumpFile[128], logFile[128];

/**
## Boundary conditions

Poiseuille inflow of the continuous phase on the left (mean velocity
`Ca` in capillary units, ramped over `tRamp`), outflow on the right,
no-slip on the embedded tube wall. The bottom boundary is the axis of
symmetry (handled by `axi.h`).
*/
#define INLET_RAMP (tRamp > 0. ? min (t/tRamp, 1.) : 1.)

u.n[left] = dirichlet (y < Rtube ?
                       2.*Ca*INLET_RAMP*(1. - sq(y/Rtube)) : 0.);
u.t[left] = dirichlet (0.);
p[left]   = neumann (0.);
pf[left]  = neumann (0.);
f[left]   = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right]   = dirichlet (0.);
pf[right]  = dirichlet (0.);

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

/**
### main()

Loads runtime parameters, sets fluid properties from the dimensionless
groups and enters the event loop.
*/
int main (int argc, char const *argv[])
{
  params_init_from_argv (argc, argv);

  CaseNo   = param_int ("CaseNo", 1000);
  MAXlevel = param_int ("MAXlevel", 10);
  MINlevel = param_int ("MINlevel", 4);

  Ca      = param_double ("Ca", 0.05);
  La      = param_double ("La", 1.);
  muR     = param_double ("muR", 1e-2);
  rhoR    = param_double ("rhoR", 1e-3);
  Rtube   = param_double ("Rtube", 0.7);
  Rb0frac = param_double ("Rb0frac", 0.8);
  xRear   = param_double ("xRear", 1.0);
  Ldomain = param_double ("Ldomain", 16.);
  /**
  The film reaches its steady thickness only after the bubble has
  travelled several of its own lengths, and the mean inlet velocity is
  $U = Ca$ in these units, so one fixed `tmax` can only ever be right at
  one capillary number. The default therefore derives from a *travel
  distance*: advancing `travelR` radii takes $t = travelR/Ca$. An
  explicit `tmax` still wins. The snapshot cadence follows it, so a case
  writes a fixed number of snapshots whatever its duration rather than
  5000 of them at low $Ca$. */

  /**
  Convergence stop. The run ends when the film thickness stops changing,
  not at an arbitrary time: `bFilm` is averaged over successive windows
  of `advWin` radii of front-tip advance, and the run stops once
  consecutive window means agree to within `bTol` for `convHold` windows,
  after at least `advMin` of advance. `tmax` remains a hard cap. */

  bTol     = param_double ("bTol", 2e-3);
  advWin   = param_double ("advWin", 0.25);
  advMin   = param_double ("advMin", 1.0);
  convHold = param_int    ("convHold", 3);

  uRel   = param_double ("uRel", 1e-2);
  dRel   = param_double ("dRel", 1e-2);
  csErr  = param_double ("csErr", 1e-2);
  VelErr = uRel*Ca;
  DErr   = dRel*Ca/Rtube;

  travelR = param_double ("travelR", 10.);
  tmax    = param_double ("tmax", travelR/Ca);
  tsnap   = param_double ("tsnap", tmax/200.);
  tRamp   = param_double ("tRamp", 1.);

  /**
  The time-step cap goes into `DT`: `centered.h` resets `dtmax = DT`
  every iteration, so assigning `dtmax` directly has no effect. */

  DT = param_double ("dtmax", 1e-2);

  /**
  The initial shape is a capsule (cylinder of length `Lcyl` with
  hemispherical caps of radius `Rb0`) whose volume equals that of the
  unit volume-equivalent sphere: $L_{cyl} = \tfrac{4}{3}(1 -
  R_{b0}^3)/R_{b0}^2$. */

  Rb0  = Rb0frac*Rtube;
  Lcyl = 4.*(1. - cube(Rb0))/(3.*sq(Rb0));
  Xb0  = xRear + Rb0 + Lcyl/2.;
  vol0 = 4.*pi/3.;

  if (CaseNo < 1000 || MAXlevel <= 0 || MINlevel <= 0 ||
      MINlevel > MAXlevel || Ca <= 0. || La <= 0. || muR <= 0. ||
      rhoR <= 0. || Rtube <= 0. || Rtube >= 1. ||
      Rb0frac <= 0. || Rb0frac >= 1. || xRear <= 0. || Ldomain <= 0. ||
      tmax <= 0. || tsnap <= 0. || DT <= 0. || tRamp < 0. ||
      travelR <= 0. || bTol <= 0. || advWin <= 0. || advMin < 0. ||
      convHold < 1 || uRel <= 0. || dRel <= 0. || csErr <= 0.) {
    fprintf (ferr, "ERROR: Invalid runtime parameters.\n");
    return 1;
  }
  if (Xb0 + Lcyl/2. + Rb0 >= Ldomain) {
    fprintf (ferr, "ERROR: initial capsule crosses the outlet; "
             "increase Ldomain or reduce xRear.\n");
    return 1;
  }
  if (Xb0 + Lcyl/2. + Rb0 > Ldomain - 4.*Rtube)
    fprintf (ferr, "WARNING: little travel room ahead of the front tip; "
             "increase Ldomain or reduce xRear.\n");

  /**
  A run asking for more advance than the domain holds would drive the
  bubble into the outlet before `tmax`. */

  double room = Ldomain - (Xb0 + Lcyl/2. + Rb0);
  if (Ca*tmax > room)
    fprintf (ferr, "WARNING: requested advance %g R exceeds the %g R of "
             "travel room ahead of the front tip; the bubble reaches the "
             "outlet before tmax = %g.\n", Ca*tmax, room, tmax);

  L0 = Ldomain;
  init_grid (1 << MINlevel);

  system ("mkdir -p intermediate");
  sprintf (dumpFile, "restart");
  sprintf (logFile, "c%d-log", CaseNo);

  /**
  Fluid 1 (`f = 1`) is the dispersed drop/bubble; fluid 2 (`f = 0`) is
  the continuous wetting phase. With $\sigma = \mu_c = R = 1$:
  $\rho_c = La$, $\rho_d = La\,\rho_R$, $\mu_d = \mu_R$. */

  rho1 = La*rhoR; mu1 = muR;
  rho2 = La;      mu2 = 1.;
  f.sigma = 1.;

  TOLERANCE = 1e-4;
  CFL = 0.5;

  if (pid() == 0) {
    fprintf (ferr, "CaseNo=%d MAXlevel=%d MINlevel=%d Ca=%g La=%g muR=%g "
             "rhoR=%g Rtube=%g Rb0=%g Lcyl=%g tmax=%g\n",
             CaseNo, MAXlevel, MINlevel, Ca, La, muR, rhoR, Rtube, Rb0,
             Lcyl, tmax);
    fprintf (ferr, "Logging to %s\n", logFile);
  }

  run();
}

/**
## Event: initialisation

Restores from `restart` when present; otherwise refines the tube
interior, embeds the wall and initialises the capsule interface. The
axi+embed metric must be resynchronised in both branches (see
`src-local/embed-vof-tube.h`).
*/
event init (t = 0)
{
  if (!restore (file = dumpFile)) {
    refine (y < 1.05*Rtube && level < MAXlevel - 2);
    refine (y < 1.05*Rtube &&
            x > Xb0 - Lcyl/2. - Rb0 - 4.*Rtube &&
            x < Xb0 + Lcyl/2. + Rb0 + 4.*Rtube && level < MAXlevel);
    tube_solid (Rtube);
    fraction (f, Rb0 - sqrt (sq (max (fabs (x - Xb0) - Lcyl/2., 0.))
                             + sq(y)));
    vof_solid_cleanup (f);
  }
  else
    embed_axi_metric_sync();
}

/**
## Adaptive mesh refinement

Adapts on the interface, its curvature, the embedded fraction (which
keeps the wall and film region refined) and the velocity, then restores
the axi+embed metric and the solid-cell volume fraction.
*/
event adapt (i++)
{
  scalar KAPPA[];
  curvature (f, KAPPA);

  /**
  The strain-rate magnitude $\sqrt{\mathbf{D}\!:\!\mathbf{D}}$ is
  adapted on directly, so the mesh follows viscous dissipation rather
  than only the interface and the velocity. In the film the shear is far
  above the bulk $U/R_{tube}$, which is exactly where the resolution is
  wanted. */

  scalar Dmag[];
  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    double D22 = (y > 1e-10) ? u.y[]/y : 0.;
    double D33 = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    double D13 = 0.5*((u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2.*Delta));
    Dmag[] = sqrt (sq(D11) + sq(D22) + sq(D33) + 2.*sq(D13));
  }

  adapt_wavelet ((scalar *){f, KAPPA, cs, u.x, u.y, Dmag},
                 (double[]){fErr, KErr, csErr, VelErr, VelErr, DErr},
                 MAXlevel, MINlevel);
  embed_axi_metric_sync();
  vof_solid_cleanup (f);
}

/**
## Event: writingFiles

Restart dump plus time-stamped snapshots in `intermediate/`.
*/
event writingFiles (t = 0; t += tsnap; t <= tmax)
{
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

/**
## Event: logWriting

Per-step diagnostics: kinetic energy, dispersed-phase volume error,
front/rear tip positions and the instantaneous film thickness
$b = R_{tube} - \max_y(\text{interface})$. The tip velocity — and hence
$Ca_b$ — is recovered in post-processing by differentiating `xTipF`.
*/
event logWriting (i++)
{
  double ke = 0.;
  foreach (reduction(+:ke))
    ke += 2.*pi*cm[]*0.5*rho(f[])*(sq(u.x[]) + sq(u.y[]))*sq(Delta);

  double vol = 2.*pi*statsf(f).sum;

  scalar xpos[], ypos[];
  position (f, xpos, {1, 0});
  position (f, ypos, {0, 1});
  double xTipF = statsf(xpos).max, xTipR = statsf(xpos).min;
  double yMax = statsf(ypos).max;
  double bFilm = Rtube - yMax;

  if (pid() == 0) {
    static FILE * fp = NULL;
    if (i == 0) {
      fp = fopen (logFile, "w");
      if (fp == NULL) {
        fprintf (ferr, "ERROR: cannot open log file %s\n", logFile);
        exit (1);
      }
      fprintf (fp, "# CaseNo %d, MAXlevel %d, Ca %g, La %g, muR %g, "
               "rhoR %g, Rtube %g\n",
               CaseNo, MAXlevel, Ca, La, muR, rhoR, Rtube);
      fprintf (fp, "# i dt t ke dVol/Vol0 xTipF xTipR bFilm\n");
    }
    else
      fp = fopen (logFile, "a");
    if (fp != NULL) {
      fprintf (fp, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
               i, dt, t, ke, (vol - vol0)/vol0, xTipF, xTipR, bFilm);
      fclose (fp);
      fp = NULL;
    }
    fprintf (ferr, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
             i, dt, t, ke, (vol - vol0)/vol0, xTipF, xTipR, bFilm);
  }

  /**
  Hard failure guards: kinetic-energy blow-up, loss of the wetting film
  (interface reaching within one fine cell of the wall) or the front
  meniscus approaching the outlet. */

  if (ke > 1e3 && i > 10) {
    if (pid() == 0)
      fprintf (ferr, "Kinetic energy blew up. Stopping.\n");
    dump (file = dumpFile);
    return 1;
  }
  if (bFilm < 4.*L0/(1 << MAXlevel) && i > 10 && pid() == 0)
    fprintf (ferr, "WARNING: film resolved by fewer than 4 cells at t=%g; "
             "increase MAXlevel.\n", t);
  if (xTipF > L0 - 2.*Rtube) {
    if (pid() == 0)
      fprintf (ferr, "Front tip reached the outlet buffer at t=%g. "
               "Stopping.\n", t);
    dump (file = dumpFile);
    return 1;
  }

  /**
  ## Film-thickness convergence

  The observable of interest is the steady film, so the run stops when
  `bFilm` stops changing rather than after a prescribed time. `bFilm` is
  accumulated over successive windows of `advWin` radii of front-tip
  advance and consecutive window means are compared; `convHold`
  successive agreements within `bTol` end the run.

  The advance origin is the *analytic* initial front position rather
  than the first sampled value, so the criterion behaves identically on
  a restarted run. Windows in advance rather than time make the test
  independent of $Ca$: a window is the same amount of interface laid
  down at every capillary number.
  */

  static double winSum = 0., winStartAdv = 0., prevMean = -1.;
  static long winN = 0;
  static int holdCount = 0;

  double adv = xTipF - (Xb0 + Lcyl/2. + Rb0);
  winSum += bFilm; winN++;

  if (adv - winStartAdv >= advWin && winN > 0) {
    double mean = winSum/winN;
    if (prevMean > 0. && adv >= advMin) {
      double drift = fabs (mean - prevMean)/prevMean;
      if (drift < bTol)
        holdCount++;
      else
        holdCount = 0;
      if (pid() == 0)
        fprintf (ferr, "# film window: adv=%.3f b=%.6e drift=%.3e "
                 "hold=%d/%d\n", adv, mean, drift, holdCount, convHold);
      if (holdCount >= convHold) {
        if (pid() == 0)
          fprintf (ferr, "Film converged: b=%.6e (b/Rtube=%.6e) after "
                   "%.3f R of advance at t=%g; drift %.3e < %.3e over "
                   "%d windows. Stopping.\n",
                   mean, mean/Rtube, adv, t, drift, bTol, convHold);
        dump (file = dumpFile);
        return 1;
      }
    }
    prevMean = mean; winSum = 0.; winN = 0; winStartAdv = adv;
  }
}

/**
## Event: stopSimulation
*/
event stopSimulation (t = tmax)
{
  if (pid() == 0)
    fprintf (ferr, "Case %d complete: Ca %g, La %g, muR %g, rhoR %g.\n",
             CaseNo, Ca, La, muR, rhoR);
  return 1;
}
