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
`Ldomain` [16], `travelR` [8], `tmax` [`travelR`/`Ca`],
`tsnap` [`tmax`/200], `tRamp` [1], `dtmax` [0.01], `bTol` [2e-3],
`advWin` [0.25], `advMin` [1.0], `convHold` [3], `uRel` [1e-2],
`dRel` [1e-2], `csErr` [1e-2], `filmBins` [128], `filmCoverage` [0.95],
`filmFlatTol` [0.05], `freshFracMin` [0.98], `speedTol` [0.02],
`shapeTol` [0.02], `filmCells` [4], `filmMinLevel` [0],
`requireShapeSteady` [true], `filmSampleDt` [derived] and `solverTol` [1e-4].
Restart controls are `restartBurnR` [0], `regridBurnR` [0] and `freshFront`
[the analytic initial front].

`tmax` and `tsnap` derive from the capillary number by default. The
bubble advances at $U = Ca$, so a fixed run time gives a different
travel distance at every $Ca$, and the film only reaches its steady
thickness after several bubble lengths of advance. `travelR` sets that
advance in units of $R$. The run stops when the front tip
reaches the outlet buffer at $L_0 - 2R_{tube}$, which leaves about
9.0 radii of usable advance, so `travelR` is capped below that. In practice the
run normally ends earlier, after the central deposited film and front speed
pass their observation-window gates (and, by default, the complete travelling
shape diagnostic also passes). `tmax` is a failure cap rather than the
expected duration.
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "params.h"
#include "embed-vof-tube.h"
#include "central-film.h"
#include <errno.h>
#include <sys/stat.h>

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
double filmCoverage, filmFlatTol, freshFracMin, speedTol, shapeTol;
double filmCells, filmSampleDt, solverTol;
double restartBurnR, regridBurnR, freshFront;
double observationStartFront = 0., observationBurnDistance = 0.;
int convHold, filmBins, filmMinLevel;
bool requireShapeSteady;
bool restoredRun = false, logNeedsHeader = false, logHasContent = false;
bool restartFilePresent = false;
int runFailed = 0;
double latestSafetyKE = 0.;

#define CENTRAL_FILM_MAX_BINS 512
CentralFilmObserver filmObserver;
CentralFilmConfig filmConfig;
CentralFilmMeasurement latestFilmMeasurement;
bool latestFilmValid = false;
double nextFilmSampleTime = 0.;

char nameOut[128], dumpFile[128], logFile[128];

/**
Reconstructs the upper VOF interface at uniformly spaced axial stations.
Each MPI rank contributes local intersections, followed by element-wise
maximum reductions.  The resolution screen uses the largest leaf intersected
by each film column, so `film/Delta_max` is a conservative lower bound on the
number of cells across that column.
*/
static bool measure_central_film (double xRearNow, double xFrontNow,
                                  double centroidNow,
                                  CentralFilmMeasurement * measurement)
{
  if (!(xFrontNow > xRearNow) || filmBins < 1)
    return false;

  double radius[CENTRAL_FILM_MAX_BINS], localRadius[CENTRAL_FILM_MAX_BINS];
  double maxDelta[CENTRAL_FILM_MAX_BINS];
  const double dx = (xFrontNow - xRearNow)/filmBins;
  for (int j = 0; j < filmBins; j++)
    radius[j] = localRadius[j] = -HUGE, maxDelta[j] = 0.;

  foreach (serial)
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord normal = interface_normal (point, f), segment[2];
      double alpha = plane_alpha (f[], normal);
      if (facets (normal, alpha, segment) == 2) {
        double x0 = x + segment[0].x*Delta;
        double y0 = y + segment[0].y*Delta;
        double x1 = x + segment[1].x*Delta;
        double y1 = y + segment[1].y*Delta;
        if (fabs(x1 - x0) > 1e-14) {
          int first = max(0, (int) ceil((min(x0, x1) - xRearNow)/dx - .5));
          int last = min(filmBins - 1,
                         (int) floor((max(x0, x1) - xRearNow)/dx - .5));
          for (int j = first; j <= last; j++) {
            double xj = xRearNow + (j + .5)*dx;
            double rj = y0 + (xj - x0)*(y1 - y0)/(x1 - x0);
            localRadius[j] = max(localRadius[j], rj);
          }
        }
      }
    }
  memcpy (radius, localRadius, filmBins*sizeof(double));
#if _MPI
  mpi_all_reduce_array (radius, MPI_DOUBLE, MPI_MAX, filmBins);
#endif

  /** Find the coarsest cell actually crossed between interface and wall. */
  foreach (serial)
    if (cs[] > 0.)
      for (int j = 0; j < filmBins; j++)
        if (radius[j] > -HUGE/2.) {
          double xj = xRearNow + (j + .5)*dx;
          if (xj >= x - Delta/2. && xj < x + Delta/2. &&
              y + Delta/2. > radius[j] && y - Delta/2. < Rtube)
            maxDelta[j] = max(maxDelta[j], Delta);
        }
#if _MPI
  mpi_all_reduce_array (maxDelta, MPI_DOUBLE, MPI_MAX, filmBins);
#endif

  double films[CENTRAL_FILM_MAX_BINS];
  int expected = 0, covered = 0, fresh = 0, count = 0;
  double minimumCells = HUGE;
  for (int j = 0; j < filmBins; j++) {
    double xi = (j + .5)/filmBins;
    if (xi <= .3 || xi >= .7)
      continue;
    expected++;
    if (radius[j] <= -HUGE/2. || radius[j] >= Rtube)
      continue;
    covered++;
    double thickness = Rtube - radius[j];
    films[count++] = thickness;
    if (xRearNow + (j + .5)*dx > freshFront)
      fresh++;
    if (maxDelta[j] > 0.)
      minimumCells = min(minimumCells, thickness/maxDelta[j]);
  }
  if (!count || !expected)
    return false;

  double sorted[CENTRAL_FILM_MAX_BINS];
  memcpy (sorted, films, count*sizeof(double));
  double median = central_film_median (sorted, count);
  memcpy (sorted, films, count*sizeof(double));
  double q05 = central_film_quantile (sorted, count, .05);
  memcpy (sorted, films, count*sizeof(double));
  double q95 = central_film_quantile (sorted, count, .95);

  *measurement = (CentralFilmMeasurement) {
    .time = t, .front = xFrontNow, .rear = xRearNow,
    .centroid = centroidNow, .film = median,
    .spatial_spread = median > 0. ? (q95 - q05)/median : HUGE,
    .coverage = (double) covered/expected,
    .fresh_fraction = covered ? (double) fresh/covered : 0.,
    .minimum_cells = minimumCells < HUGE ? minimumCells : 0.
  };
  return true;
}

/**
Applies two independent controls. `filmMinLevel` floors the complete bubble
interface and its wall-facing liquid band, so changing it is a genuine grid
study of deposition and both menisci. `filmCells` may demand a still finer
level in the measured central film. Neither control refines unrelated parts of
the square domain or the solid exterior.
*/
static void refine_measured_film (void)
{
  if (!latestFilmValid || (filmCells <= 0. && filmMinLevel <= 0) ||
      !(latestFilmMeasurement.film > 0.))
    return;
  int measuredTarget = 0;
  if (filmCells > 0.)
    measuredTarget = (int) ceil
      (log(L0*filmCells/latestFilmMeasurement.film)/log(2.));
  measuredTarget = max(MINlevel, min(MAXlevel, measuredTarget));
  int floorLevel = max(MINlevel, min(MAXlevel, filmMinLevel));
  const double length = latestFilmMeasurement.front - latestFilmMeasurement.rear;
  const double central0 = latestFilmMeasurement.rear + .25*length;
  const double central1 = latestFilmMeasurement.rear + .75*length;
  const double wall0 = latestFilmMeasurement.rear - Rtube;
  const double wall1 = latestFilmMeasurement.front + Rtube;
  const double band = 2.*latestFilmMeasurement.film;
  tube_refinement_geometry_begin();
  refine ((filmMinLevel > 0 && level < floorLevel && cs[] > 0. &&
           f[] > 1e-6 && f[] < 1. - 1e-6) ||
          (filmMinLevel > 0 && level < floorLevel && cs[] > 0. &&
           x > wall0 && x < wall1 && y + Delta/2. > Rtube - band) ||
          (filmCells > 0. && level < measuredTarget && cs[] > 0. &&
           x > central0 && x < central1 &&
           y + Delta/2. > Rtube - 1.5*latestFilmMeasurement.film));
  tube_refinement_geometry_end();
}

/** Returns interface tips, global minimum gap and dispersed-phase centroid. */
static bool current_drop_geometry (double * xFrontNow, double * xRearNow,
                                   double * minimumFilm,
                                   double * centroidNow)
{
  scalar xpos[], ypos[];
  position (f, xpos, {1, 0});
  position (f, ypos, {0, 1});
  *xFrontNow = statsf(xpos).max;
  *xRearNow = statsf(xpos).min;
  *minimumFilm = Rtube - statsf(ypos).max;

  double volume = 0., moment = 0.;
  foreach (reduction(+:volume) reduction(+:moment)) {
    double element = f[]*dv();
    volume += element;
    moment += x*element;
  }
  *centroidNow = volume > 0. ? moment/volume : 0.;
  return *xFrontNow > *xRearNow && volume > 0.;
}

static bool close_parameter (double a, double b)
{
  return fabs(a - b) <= 1e-12*max(1., max(fabs(a), fabs(b)));
}

/**
An existing log may only be continued when its case metadata and central-film
schema match. Missing or empty logs receive exactly one header.
*/
static bool prepare_case_log (void)
{
  struct stat status;
  if (stat(logFile, &status) != 0) {
    if (errno == ENOENT) {
      logNeedsHeader = true;
      return true;
    }
    return false;
  }
  if (status.st_size == 0) {
    logNeedsHeader = true;
    return true;
  }

  FILE * fp = fopen (logFile, "r");
  if (!fp)
    return false;
  char metadata[512] = "", columns[512] = "";
  bool readOK = fgets(metadata, sizeof(metadata), fp) &&
                fgets(columns, sizeof(columns), fp);
  fclose (fp);
  int caseNumber = 0, maximumLevel = 0, minimumLevel = 0;
  double ca = 0., la = 0., viscosityRatio = 0., densityRatio = 0., radius = 0.;
  double domain = 0., tolerance = 0.;
  bool metadataOK = readOK &&
    sscanf(metadata, "# CaseNo %d, MAXlevel %d, MINlevel %d, Ca %lf, "
           "La %lf, muR %lf, rhoR %lf, Rtube %lf, Ldomain %lf, "
           "solverTol %lf", &caseNumber, &maximumLevel, &minimumLevel, &ca,
           &la, &viscosityRatio, &densityRatio, &radius, &domain,
           &tolerance) == 10 &&
    caseNumber == CaseNo && maximumLevel == MAXlevel &&
    minimumLevel == MINlevel &&
    close_parameter(ca, Ca) && close_parameter(la, La) &&
    close_parameter(viscosityRatio, muR) &&
    close_parameter(densityRatio, rhoR) && close_parameter(radius, Rtube) &&
    close_parameter(domain, Ldomain) && close_parameter(tolerance, solverTol) &&
    strstr(columns, "bCentral") && strstr(columns, "minCells");
  if (!metadataOK)
    return false;
  logNeedsHeader = false;
  logHasContent = true;
  return true;
}

/**
An init event which returns after `restore()` leaves Basilisk to reschedule
events with the restored clock and an unset next timestep. Terminate directly
instead. All init predicates below are collective; the barrier lets rank zero
flush the exact terminal reason before MPI aborts. The existing checkpoint is
never rewritten on this path.
*/
static void terminate_init_failure (void)
{
  fflush (NULL);
#if _MPI
  MPI_Abort (MPI_COMM_WORLD, 1);
#endif
  exit (EXIT_FAILURE);
}

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
  Convergence stop. Uniform central-film measurements are averaged over
  successive windows of `advWin` radii of front-tip advance. Film and front
  speed must remain within tolerance for `convHold` windows after at least
  `advMin` of newly observed advance. `tmax` remains a hard cap. */

  bTol     = param_double ("bTol", 2e-3);
  advWin   = param_double ("advWin", 0.25);
  advMin   = param_double ("advMin", 1.0);
  convHold = param_int    ("convHold", 3);

  filmBins      = param_int    ("filmBins", 128);
  filmCoverage  = param_double ("filmCoverage", .95);
  filmFlatTol   = param_double ("filmFlatTol", .05);
  freshFracMin  = param_double ("freshFracMin", .98);
  speedTol      = param_double ("speedTol", .02);
  shapeTol      = param_double ("shapeTol", .02);
  filmCells     = param_double ("filmCells", 4.);
  filmMinLevel  = param_int    ("filmMinLevel", 0);
  requireShapeSteady = param_bool ("requireShapeSteady", true);
  restartBurnR = param_double ("restartBurnR", 0.);
  regridBurnR  = param_double ("regridBurnR", 0.);

  uRel   = param_double ("uRel", 1e-2);
  dRel   = param_double ("dRel", 1e-2);
  csErr  = param_double ("csErr", 1e-2);
  VelErr = uRel*Ca;
  DErr   = dRel*Ca/Rtube;

  travelR = param_double ("travelR", 8.);
  tmax    = param_double ("tmax", travelR/Ca);
  tsnap   = param_double ("tsnap", tmax/200.);
  tRamp   = param_double ("tRamp", 1.);
  filmSampleDt = param_double
    ("filmSampleDt", min(tsnap, advWin/(8.*Ca)));

  /**
  The time-step cap goes into `DT`: `centered.h` resets `dtmax = DT`
  every iteration, so assigning `dtmax` directly has no effect. */

  DT = param_double ("dtmax", 1e-2);
  solverTol = param_double ("solverTol", 1e-4);

  /**
  The initial shape is a capsule (cylinder of length `Lcyl` with
  hemispherical caps of radius `Rb0`) whose volume equals that of the
  unit volume-equivalent sphere: $L_{cyl} = \tfrac{4}{3}(1 -
  R_{b0}^3)/R_{b0}^2$. */

  Rb0  = Rb0frac*Rtube;
  Lcyl = 4.*(1. - cube(Rb0))/(3.*sq(Rb0));
  Xb0  = xRear + Rb0 + Lcyl/2.;
  vol0 = 4.*pi/3.;
  freshFront = param_double ("freshFront", Xb0 + Lcyl/2. + Rb0);

  const double finiteParameters[] = {
    Ca, La, muR, rhoR, Rtube, Rb0frac, xRear, Ldomain, travelR,
    tmax, tsnap, DT, tRamp, bTol, advWin, advMin, uRel, dRel, csErr,
    filmCoverage, filmFlatTol, freshFracMin, speedTol, shapeTol,
    filmCells, filmSampleDt, solverTol, restartBurnR, regridBurnR, freshFront
  };
  for (unsigned int j = 0; j < sizeof(finiteParameters)/sizeof(double); j++)
    if (!isfinite(finiteParameters[j])) {
      fprintf (ferr, "ERROR: Runtime parameters must be finite.\n");
      return 1;
    }

  if (CaseNo < 1000 || MAXlevel <= 0 || MINlevel <= 0 ||
      MINlevel > MAXlevel || Ca <= 0. || La <= 0. || muR <= 0. ||
      rhoR <= 0. || Rtube <= 0. || Rtube >= 1. ||
      Rb0frac <= 0. || Rb0frac >= 1. || xRear <= 0. || Ldomain <= 0. ||
      tmax <= 0. || tsnap <= 0. || DT <= 0. || tRamp < 0. ||
      travelR <= 0. || bTol <= 0. || advWin <= 0. || advMin < 0. ||
      convHold < 1 || uRel <= 0. || dRel <= 0. || csErr <= 0. ||
      filmBins < 32 || filmBins > CENTRAL_FILM_MAX_BINS ||
      filmCoverage <= 0. || filmCoverage > 1. || filmFlatTol <= 0. ||
      freshFracMin < 0. || freshFracMin > 1. || speedTol <= 0. ||
      shapeTol <= 0. || filmCells < 0. || filmMinLevel < 0 ||
      filmMinLevel > MAXlevel || filmSampleDt <= 0. || solverTol <= 0. ||
      restartBurnR < 0. || regridBurnR < 0. || !isfinite(freshFront)) {
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

  double room = (Ldomain - 2.*Rtube) - (Xb0 + Lcyl/2. + Rb0);
  if (Ca*tmax > room)
    fprintf (ferr, "WARNING: requested advance %g R exceeds the %g R of "
             "travel room ahead of the front tip; the bubble reaches the "
             "outlet before tmax = %g.\n", Ca*tmax, room, tmax);

  L0 = Ldomain;
  init_grid (1 << MINlevel);

  int directoryError = 0;
  if (pid() == 0 && mkdir ("intermediate", 0777) != 0) {
    int mkdirError = errno;
    struct stat status;
    if (mkdirError != EEXIST || stat ("intermediate", &status) != 0 ||
        !S_ISDIR(status.st_mode)) {
      errno = mkdirError;
      perror ("intermediate");
      directoryError = 1;
    }
  }
#if _MPI
  mpi_all_reduce (directoryError, MPI_INT, MPI_MAX);
#endif
  if (directoryError)
    return 1;
  sprintf (dumpFile, "restart");
  sprintf (logFile, "c%d-log", CaseNo);

  int restartPresent = 0;
  if (pid() == 0) {
    struct stat restartStatus;
    restartPresent = stat(dumpFile, &restartStatus) == 0 &&
                     restartStatus.st_size > 0;
  }

  int logSetupError = 0;
  if (pid() == 0 && !prepare_case_log()) {
    fprintf (ferr, "HARDFAIL_LOG_SCHEMA: existing log %s is incompatible with this case "
             "or cannot be read; refusing to append a mixed case log.\n",
             logFile);
    logSetupError = 1;
  }
  int headerNeeded = logNeedsHeader, existingLog = logHasContent;
#if _MPI
  mpi_all_reduce (restartPresent, MPI_INT, MPI_MAX);
  mpi_all_reduce (logSetupError, MPI_INT, MPI_MAX);
  mpi_all_reduce (headerNeeded, MPI_INT, MPI_MAX);
  mpi_all_reduce (existingLog, MPI_INT, MPI_MAX);
#endif
  restartFilePresent = restartPresent;
  logNeedsHeader = headerNeeded;
  logHasContent = existingLog;
  if (logSetupError)
    return 1;

  /**
  Fluid 1 (`f = 1`) is the dispersed drop/bubble; fluid 2 (`f = 0`) is
  the continuous wetting phase. With $\sigma = \mu_c = R = 1$:
  $\rho_c = La$, $\rho_d = La\,\rho_R$, $\mu_d = \mu_R$. */

  rho1 = La*rhoR; mu1 = muR;
  rho2 = La;      mu2 = 1.;
  f.sigma = 1.;

  filmConfig = (CentralFilmConfig) {
    .advance_window = advWin,
    .minimum_observed_advance = advMin,
    .film_relative_tolerance = bTol,
    .front_speed_relative_tolerance = speedTol,
    .minimum_coverage = filmCoverage,
    .maximum_spatial_spread = filmFlatTol,
    .minimum_fresh_fraction = freshFracMin,
    .minimum_cells = filmCells,
    .shape_relative_tolerance = shapeTol,
    .hold_windows = convHold
  };
  central_film_observer_init (&filmObserver, filmConfig);

  TOLERANCE = solverTol;
  CFL = 0.5;

  if (pid() == 0) {
    fprintf (ferr, "CaseNo=%d MAXlevel=%d MINlevel=%d Ca=%g La=%g muR=%g "
             "rhoR=%g Rtube=%g Rb0=%g Lcyl=%g tmax=%g solverTol=%g\n",
             CaseNo, MAXlevel, MINlevel, Ca, La, muR, rhoR, Rtube, Rb0,
             Lcyl, tmax, solverTol);
    fprintf (ferr, "Central film: bins=%d sampleDt=%g coverage>=%g "
             "spread<=%g fresh>=%g cells>=%g filmMinLevel=%d "
             "requireShapeSteady=%d restartBurnR=%g regridBurnR=%g "
             "freshFront=%g\n", filmBins, filmSampleDt,
             filmCoverage, filmFlatTol, freshFracMin, filmCells,
             filmMinLevel, requireShapeSteady, restartBurnR, regridBurnR,
             freshFront);
    fprintf (ferr, "Logging to %s\n", logFile);
  }

  run();
  return runFailed ? 1 : 0;
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
  restoredRun = restore (file = dumpFile);
  if (!restoredRun && restartFilePresent) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_RESTORE: restart exists but could not be "
               "restored; refusing to replace it with a fresh case.\n");
    terminate_init_failure();
  }
  if (!restoredRun && logHasContent) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_LOG_REUSE: %s contains an existing case log "
               "but no compatible restart was restored.\n", logFile);
    terminate_init_failure();
  }
  if (!restoredRun) {
    refine (y < 1.05*Rtube && level < MAXlevel - 2);
    refine (y < 1.05*Rtube &&
            x > Xb0 - Lcyl/2. - Rb0 - 4.*Rtube &&
            x < Xb0 + Lcyl/2. + Rb0 + 4.*Rtube && level < MAXlevel);
    tube_solid (Rtube);
    fraction (f, Rb0 - sqrt (sq (max (fabs (x - Xb0) - Lcyl/2., 0.))
                             + sq(y)));
    vof_solid_cleanup (f);
  }
  else {
    if (!close_parameter(L0, Ldomain)) {
      if (pid() == 0)
        fprintf (ferr, "HARDFAIL_DOMAIN_MISMATCH: restart L0=%g differs "
                 "from requested Ldomain=%g; a dump cannot be extended or "
                 "shrunk in place.\n", L0, Ldomain);
      terminate_init_failure();
    }
    // Dumps contain cell fields, not the embedded face fractions. Rebuild
    // the stationary wall before the solver initializes face fluxes.
    tube_solid (Rtube);
#if TREE && AXI
    // VOF's face scan also visits inactive full-fluid parent cells.
    // Initialize their exact axisymmetric metric, not only active leaves.
    foreach_cell()
      if (cs[] >= 1.)
        cm[] = y;
#endif
    vof_solid_cleanup (f);
  }

  /**
  A warm restart may request a finer film floor or a measured cell target.
  Reconstruct the restored geometry, then perform one bounded recursive
  refinement to the requested level (never beyond `MAXlevel`).  Observation
  history is intentionally not restored: the run must cover new advance
  windows before any terminal success.
  */
  double front, rear, minimumFilm, centroid;
  if (!current_drop_geometry (&front, &rear, &minimumFilm, &centroid)) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_GEOMETRY: cannot reconstruct the dispersed "
               "phase after initialisation.\n");
    terminate_init_failure();
  }
  observationStartFront = front;
  observationBurnDistance = restoredRun ? max(restartBurnR, regridBurnR) : 0.;
  central_film_observer_reset (&filmObserver);
  bool initialFilmOK = measure_central_film
    (rear, front, centroid, &latestFilmMeasurement);
  if (initialFilmOK)
    latestFilmValid = true;
  else if (filmMinLevel > 0 || filmCells > 0.) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_REGRID_SAMPLE: cannot measure the film needed "
               "to apply filmMinLevel/filmCells before integration.\n");
    terminate_init_failure();
  }
  refine_measured_film();
  embed_axi_metric_sync();
  vof_solid_cleanup (f);
  if (pid() == 0 && restoredRun)
    fprintf (ferr, "# restart observation burn: startFront=%g distance=%g "
             "freshFront=%g\n", observationStartFront,
             observationBurnDistance, freshFront);
}

/**
Cheap per-step corruption guard. Interface reconstruction stays on the
bounded diagnostic cadence, while nonfinite primary fields, time-step failure
and kinetic-energy blow-up terminate before several corrupt steps accumulate.
*/
event solverSafety (i++)
{
  double ke = 0.;
  double nonfinite = (!isfinite(t) || !isfinite(dt) || dt <= 0.) ? 1. : 0.;
  foreach (reduction(+:ke) reduction(+:nonfinite)) {
    if (!isfinite(f[]) || !isfinite(u.x[]) || !isfinite(u.y[]) ||
        !isfinite(p[]) || !isfinite(cs[]) || !isfinite(cm[]))
      nonfinite += 1.;
    ke += 2.*pi*cm[]*0.5*rho(f[])*(sq(u.x[]) + sq(u.y[]))*sq(Delta);
  }
  latestSafetyKE = ke;
  if (nonfinite || !isfinite(ke)) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_NONFINITE: nonfinite solver state at i=%d "
               "t=%g dt=%g.\n", i, t, dt);
    runFailed = 1;
    if (pid() == 0)
      fprintf (ferr, "# last good restart preserved at %s.\n", dumpFile);
    return 1;
  }
  if (ke > 1e3 && i > 10) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_ENERGY: kinetic energy %.6e exceeds 1e3 "
               "at i=%d t=%g.\n", ke, i, t);
    runFailed = 1;
    dump (file = "hardfail-energy");
    return 1;
  }
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
    double D22 = (y > 1e-10) ? u.y[]/y : D11;  // axis limit is du_r/dr
    double D33 = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    double D13 = 0.5*((u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2.*Delta));
    Dmag[] = sqrt (sq(D11) + sq(D22) + sq(D33) + 2.*sq(D13));
  }

  adapt_wavelet ((scalar *){f, KAPPA, cs, u.x, u.y, Dmag},
                 (double[]){fErr, KErr, csErr, VelErr, VelErr, DErr},
                 MAXlevel, MINlevel);
  refine_measured_film();
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

Bounded-cadence diagnostics: kinetic energy, dispersed-phase volume error,
front/rear/centroid positions, the legacy global minimum gap
$b_{min} = R_{tube} - \max_y(\text{interface})$, and the uniformly sampled
central film. The window diagnostics directly report front, rear and centroid
speeds.
*/
event logWriting (i++)
{
  /** Interface reconstruction is deliberately sampled at a bounded cadence. */
  if (t + 1e-12 < nextFilmSampleTime)
    return 0;
  nextFilmSampleTime = t + filmSampleDt;

  double ke = latestSafetyKE;
  double vol = 2.*pi*statsf(f).sum;
  double xTipF, xTipR, bFilm, centroid;
  bool geometryOK = current_drop_geometry
    (&xTipF, &xTipR, &bFilm, &centroid);
  bool centralOK = geometryOK && measure_central_film
    (xTipR, xTipF, centroid, &latestFilmMeasurement);
  if (centralOK)
    latestFilmValid = true;
  else
    central_film_observer_reset (&filmObserver);

  double central = centralOK ? latestFilmMeasurement.film : NAN;
  double spread = centralOK ? latestFilmMeasurement.spatial_spread : NAN;
  double coverage = centralOK ? latestFilmMeasurement.coverage : 0.;
  double fresh = centralOK ? latestFilmMeasurement.fresh_fraction : 0.;
  double cells = centralOK ? latestFilmMeasurement.minimum_cells : 0.;

  int logError = 0;
  if (pid() == 0) {
    FILE * fp = fopen (logFile, "a");
    if (fp == NULL) {
      fprintf (ferr, "HARDFAIL_LOG_IO: cannot open log file %s\n", logFile);
      logError = 1;
    }
    else {
      if (logNeedsHeader) {
        fprintf (fp, "# CaseNo %d, MAXlevel %d, MINlevel %d, Ca %.17g, "
                 "La %.17g, muR %.17g, rhoR %.17g, Rtube %.17g, "
                 "Ldomain %.17g, solverTol %.17g\n", CaseNo, MAXlevel,
                 MINlevel, Ca, La, muR, rhoR, Rtube, Ldomain, solverTol);
        fprintf (fp, "# i dt t ke dVol/Vol0 xTipF xTipR bFilm "
                 "bCentral spatialSpread coverage freshFraction minCells "
                 "centroid length\n");
        logNeedsHeader = false;
      }
      fprintf (fp, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e "
               "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
               i, dt, t, ke, (vol - vol0)/vol0, xTipF, xTipR, bFilm,
               central, spread, coverage, fresh, cells, centroid,
               xTipF - xTipR);
      fclose (fp);
    }
    fprintf (ferr, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e "
             "%.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
             i, dt, t, ke, (vol - vol0)/vol0, xTipF, xTipR, bFilm,
             central, spread, coverage, fresh, cells, centroid,
             xTipF - xTipR);
  }
  if (logError) {
#if _MPI
    // A rank-zero I/O failure must end this MPI job, without adding a
    // collective reduction to every successful diagnostic step.
    MPI_Abort (MPI_COMM_WORLD, 1);
#endif
    runFailed = 1;
    return 1;
  }

  /**
  Hard failure guards: kinetic-energy blow-up, loss of the wetting film
  (interface reaching within one fine cell of the wall) or the front
  meniscus approaching the outlet. */

  if (bFilm < 4.*L0/(1 << MAXlevel) && i > 10 && pid() == 0)
    fprintf (ferr, "WARNING: film resolved by fewer than 4 cells at t=%g; "
             "increase MAXlevel.\n", t);
  if (xTipF > L0 - 2.*Rtube) {
    if (pid() == 0)
      fprintf (ferr, "INCOMPLETE_OUTLET: front tip reached the outlet "
               "buffer at t=%g.\n", t);
    runFailed = 1;
    dump (file = dumpFile);
    return 1;
  }

  /**
  Central-film stationarity is based on uniformly sampled 30--70% profiles,
  successive advance windows and a separately measured front-speed plateau.
  Whole-shape stationarity compares front, rear and centroid speeds plus
  length drift; it is a distinct milestone and can optionally be required.
  */
  bool observationEligible = centralOK && t >= tRamp &&
    xTipF >= observationStartFront + observationBurnDistance;
  if (!observationEligible)
    central_film_observer_reset (&filmObserver);
  if (observationEligible) {
    CentralFilmWindow window;
    if (central_film_observer_push
        (&filmObserver, latestFilmMeasurement, &window)) {
      if (pid() == 0)
        fprintf (ferr, "# central-film window: advance=%.3f b=%.6e "
                 "filmDrift=%.3e frontSpeed=%.6e speedDrift=%.3e "
                 "coverage=%.3f spread=%.3e fresh=%.3f minCells=%.3f "
                 "shapeMismatch=%.3e lengthDrift=%.3e filmHold=%d/%d "
                 "shapeHold=%d/%d\n",
                 window.advance, window.film, window.film_drift,
                 window.front_speed, window.front_speed_drift,
                 window.minimum_coverage, window.maximum_spatial_spread,
                 window.minimum_fresh_fraction, window.minimum_cells,
                 window.shape_speed_mismatch, window.length_drift,
                 window.hold_count, convHold, window.shape_hold_count,
                 convHold);

      if (filmCells > 0. && window.minimum_cells < filmCells && pid() == 0)
        fprintf (ferr, "# central-film resolution insufficient: measured "
                 "minimum %.3f cells, requested %.3f; refinement is capped "
                 "at MAXlevel=%d.\n", window.minimum_cells, filmCells,
                 MAXlevel);

      if (window.film_converged) {
        if (pid() == 0)
          fprintf (ferr, "Central deposited-film steady: b=%.6e "
                   "(b/Rtube=%.6e), front speed %.6e. Whole-shape status: "
                   "%s (%d/%d consecutive windows).\n", window.film,
                   window.film/Rtube,
                   window.front_speed,
                   window.shape_converged ? "steady" : "still evolving",
                   window.shape_hold_count, convHold);
        if (!requireShapeSteady || window.shape_converged) {
          if (pid() == 0)
            fprintf (ferr, "SUCCESS: %s milestone passed at t=%g.\n",
                     requireShapeSteady ? "central-film and whole-shape" :
                                          "central deposited-film", t);
          dump (file = dumpFile);
          return 1;
        }
      }
    }
  }
}

/**
## Event: stopSimulation
*/
event stopSimulation (t = tmax)
{
  runFailed = 1;
  if (pid() == 0)
    fprintf (ferr, "INCOMPLETE_TMAX: case %d reached tmax without the requested "
             "stationarity milestone: Ca %g, La %g, muR %g, rhoR %g.\n",
             CaseNo, Ca, La, muR, rhoR);
  dump (file = dumpFile);
  return 1;
}
