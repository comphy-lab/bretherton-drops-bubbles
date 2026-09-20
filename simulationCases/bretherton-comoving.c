/**
# bretherton-comoving.c

Axisymmetric two-phase simulation of a long bubble in a liquid-filled
capillary tube, computed in a window that translates with the bubble.
The physics, non-dimensionalisation, embedded wall and VOF machinery are
those of [`bretherton.c`](bretherton.c); this case changes the frame.

## Frame

Lab velocity $\mathbf{v} = \mathbf{u} + U(t)\mathbf{e}_x$. The window
equations are the same Navier--Stokes equations plus the uniform body force
$-\dot U\,\mathbf{e}_x$ per unit mass on both phases (an exact Galilean
transformation). The tube wall moves at $-U$, and both ends carry the fully
developed Poiseuille profile shifted by $-U$, which is the infinite-tube
condition in a circular tube. $U$ is set by a critically damped controller
on the dispersed-phase centroid (`src-local/comoving-frame.h`); Galilean
invariance decouples it from the physics, so $U \to U_b = Ca_b$ at steady
state and the reported bubble speed is a measured frame speed, not a fit
to a tip trajectory.

In the frame the mean liquid speed is $Ca_{in} - U_b < 0$, so liquid enters
through the *front* end and leaves through the *rear* end. With a velocity
Dirichlet condition on every boundary the pressure is defined up to a
constant; pressures are reported as gauges relative to the mean pressure on
the front end. The discrete inflow and outflow match only when the leaf
columns at both ends share the same level, so both end columns are held at
`padLevel`. The residual mismatch is logged.

## Steady state

Judged in film renewals $L_b/U$ (`wRenew` window, `burnRenew` burn-in,
`convHold` consecutive stable windows). The film freshness gate reuses the
parent's `freshFront` with the frame displacement subtracted.

## Runtime keys

As `bretherton.c` plus: `tau` [1], `Uframe0` [0], `xTarget` [initial
centroid], `prescribedU` [false], `padLevel` [MAXlevel-2], `wRenew` [0.5],
`burnRenew` [1.0], `posTol` [0.02] (in `Rtube`), `CaPrev` [0], `outlet`
[`poiseuille`|`open-rear`|`parent`].
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
#include "comoving-frame.h"
#include "tag.h"
#include <errno.h>
#include <sys/stat.h>

#define fErr (1e-3)
#define KErr (1e-4)

int MAXlevel, MINlevel, CaseNo, padLevel;
double Ca, CaPrev, La, muR, rhoR, Rtube, Rb0frac, xRear, Ldomain;
double tmax = 200., tsnap = 1., tRamp = 1., travelR = 10.;
double Rb0, Lcyl, Xb0, vol0;
double bTol, uRel, dRel, VelErr, DErr, csErr;
double filmCoverage, filmFlatTol, freshFracMin, speedTol, shapeTol, posTol;
double filmCells, filmSampleDt, solverTol;
double freshFront, tau, Uframe0, xTarget, wRenew, burnRenew;
int convHold, filmBins, filmMinLevel;
bool prescribedU;
bool restoredRun = false, logNeedsHeader = false, logHasContent = false;
bool restartFilePresent = false;
int runFailed = 0, outletMode = 0;   // 0 poiseuille, 1 open-rear, 2 parent
double latestSafetyKE = 0., tStart = 0.;

#define CENTRAL_FILM_MAX_BINS 512
CentralFilmMeasurement latestFilmMeasurement;
bool latestFilmValid = false;
double nextFilmSampleTime = 0.;
ComovingFrame frame;
RenewalObserver renewal;
face vector av[];

char nameOut[128], dumpFile[128], logFile[128], frameFile[128];

/** Inlet capillary number ramped from `CaPrev` to `Ca` over `tRamp`. */
static inline double ca_now (void)
{
  double r = tRamp > 0. ? min ((t - tStart)/tRamp, 1.) : 1.;
  return CaPrev + (Ca - CaPrev)*r;
}
/**
The end profiles are imposed as the fluid-weighted mean of the parabola over
each boundary face, $\bar u = \int 2Ca(1-r^2/R_t^2)\,r\,dr / \int r\,dr$
over the fluid part of the face. With the tube face metric the discrete
boundary flux is then exact at any level, and the wall-adjacent cut cell,
whose centre lies outside the tube, still receives the $-U$ shift. A point
value sampled at the cell centre gave that cut cell zero instead of $-U$,
which in the moving frame is equivalent to a 5 % larger inlet flux. */
static inline double poiseuille_face_mean (double yc, double d)
{
  double lo = yc - d/2., hi = min (yc + d/2., Rtube);
  if (hi <= lo)
    return 0.;
  double a2 = hi*hi - lo*lo, a4 = sq(hi*hi) - sq(lo*lo);
  return 2.*ca_now()*(a2 - a4/(2.*sq(Rtube)))/a2;
}
#define POISEUILLE_FRAME (y - Delta/2. < Rtube ?                        \
                          poiseuille_face_mean (y, Delta) - frame.U : 0.)

/**
Reconstructs the upper VOF interface at uniformly spaced axial stations
(same estimator as the parent). */
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
        double x0 = x + segment[0].x*Delta, y0 = y + segment[0].y*Delta;
        double x1 = x + segment[1].x*Delta, y1 = y + segment[1].y*Delta;
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

/** Film floor and measured-film refinement, as in the parent, plus the fixed
end columns that make the two boundary fluxes discretely identical. */
static void refine_measured_film (void)
{
  const double endBand = 2.*L0/(1 << padLevel);
  tube_refinement_geometry_begin();
  refine (level < padLevel && cs[] > 0. &&
          (x < endBand || x > L0 - endBand));
  if (latestFilmValid && (filmCells > 0. || filmMinLevel > 0) &&
      latestFilmMeasurement.film > 0.) {
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
    refine ((filmMinLevel > 0 && level < floorLevel && cs[] > 0. &&
             f[] > 1e-6 && f[] < 1. - 1e-6) ||
            (filmMinLevel > 0 && level < floorLevel && cs[] > 0. &&
             x > wall0 && x < wall1 && y + Delta/2. > Rtube - band) ||
            (filmCells > 0. && level < measuredTarget && cs[] > 0. &&
             x > central0 && x < central1 &&
             y + Delta/2. > Rtube - 1.5*latestFilmMeasurement.film));
  }
  tube_refinement_geometry_end();
}

/** Cheap centroid of all dispersed phase, for the controller. */
static double bubble_centroid (void)
{
  double volume = 0., moment = 0.;
  foreach (reduction(+:volume) reduction(+:moment)) {
    double element = f[]*dv();
    volume += element;
    moment += x*element;
  }
  return volume > 0. ? moment/volume : 0.;
}

/** Tips, minimum gap, centroid and rear-cap indentation on the largest
connected VOF component. */
static bool current_drop_geometry (double * xFrontNow, double * xRearNow,
                                   double * minimumFilm, double * centroidNow,
                                   double * rearAxisTip, int * components)
{
  scalar cc[];
  foreach()
    cc[] = f[] > 1e-4;
  int n = tag (cc);
  *components = n;
  if (n < 1)
    return false;
  double volComp[n];
  for (int j = 0; j < n; j++)
    volComp[j] = 0.;
  foreach (serial)
    if (cc[] > 0.)
      volComp[((int) cc[]) - 1] += f[]*dv();
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, volComp, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  int mainTag = 1;
  double vMain = volComp[0];
  for (int j = 1; j < n; j++)
    if (volComp[j] > vMain) {
      vMain = volComp[j];
      mainTag = j + 1;
    }
  scalar xpos[], ypos[];
  position (f, xpos, {1, 0});
  position (f, ypos, {0, 1});
  double xFront = -HUGE, xRear = HUGE, yMax = -HUGE, volume = 0., moment = 0.;
  double xAxis = HUGE;
  foreach (reduction(max:xFront) reduction(min:xRear) reduction(max:yMax)
           reduction(+:volume) reduction(+:moment) reduction(min:xAxis)) {
    if ((int) cc[] != mainTag)
      continue;
    double element = f[]*dv();
    volume += element;
    moment += x*element;
    if (cs[] > 0. && f[] > 1e-6 && f[] < 1. - 1e-6) {
      if (finite (xpos[]) && xpos[] > xFront) xFront = xpos[];
      if (finite (xpos[]) && xpos[] < xRear) xRear = xpos[];
      if (finite (ypos[]) && ypos[] > yMax) yMax = ypos[];
      if (finite (xpos[]) && y < 1.5*Delta && xpos[] < xAxis) xAxis = xpos[];
    }
  }
  *xFrontNow = xFront;
  *xRearNow = xRear;
  *minimumFilm = Rtube - yMax;
  *centroidNow = volume > 0. ? moment/volume : 0.;
  *rearAxisTip = xAxis < HUGE ? xAxis : xRear;
  return finite (xFront) && finite (xRear) && xFront > xRear && volume > 0.;
}

/** Boundary fluxes (2 pi omitted) and mean front-end pressure. */
static void boundary_fluxes (double * Qrear, double * Qfront, double * pFront)
{
  double ql = 0., qr = 0., ps = 0., pw = 0.;
  foreach_boundary (left, reduction(+:ql))
    ql += uf.x[]*Delta;
  foreach_boundary (right, reduction(+:qr) reduction(+:ps) reduction(+:pw)) {
    qr += uf.x[1]*Delta;
    if (cs[] > 0.) {
      ps += p[]*cm[]*Delta;
      pw += cm[]*Delta;
    }
  }
  *Qrear = ql; *Qfront = qr;
  *pFront = pw > 0. ? ps/pw : 0.;
}

/** Mean pressure inside the bubble away from the interface. */
static double bubble_pressure (void)
{
  double ps = 0., pw = 0.;
  foreach (reduction(+:ps) reduction(+:pw))
    if (f[] > 0.99 && cs[] > 0.) {
      ps += p[]*dv();
      pw += dv();
    }
  return pw > 0. ? ps/pw : NAN;
}

static bool close_parameter (double a, double b)
{
  return fabs(a - b) <= 1e-12*max(1., max(fabs(a), fabs(b)));
}

static bool prepare_case_log (void)
{
  struct stat status;
  if (stat(logFile, &status) != 0) {
    if (errno == ENOENT) { logNeedsHeader = true; return true; }
    return false;
  }
  if (status.st_size == 0) { logNeedsHeader = true; return true; }
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
    strstr(columns, "Uframe") && strstr(columns, "Qfront");
  if (!metadataOK)
    return false;
  logNeedsHeader = false;
  logHasContent = true;
  return true;
}

static void terminate_init_failure (void)
{
  fflush (NULL);
#if _MPI
  MPI_Abort (MPI_COMM_WORLD, 1);
#endif
  exit (EXIT_FAILURE);
}

/** Frame scalars are not part of a dump; they travel in `frame-state`. */
static void write_frame_state (void)
{
  if (pid() != 0)
    return;
  FILE * fp = fopen (frameFile, "w");
  if (!fp)
    return;
  fprintf (fp, "%.17g %.17g %.17g %.17g %.17g %.17g\n", t, frame.U, frame.dUdt,
           frame.target, frame.displacement, freshFront);
  fclose (fp);
}

static bool read_frame_state (double * tt, double * U, double * dUdt,
                              double * target, double * disp, double * fresh)
{
  double v[6];
  int ok = 0;
  if (pid() == 0) {
    FILE * fp = fopen (frameFile, "r");
    if (fp) {
      ok = fscanf (fp, "%lf %lf %lf %lf %lf %lf", &v[0], &v[1], &v[2], &v[3],
                   &v[4], &v[5]) == 6;
      fclose (fp);
    }
  }
#if _MPI
  mpi_all_reduce (ok, MPI_INT, MPI_MAX);
  MPI_Bcast (v, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  if (!ok)
    return false;
  *tt = v[0]; *U = v[1]; *dUdt = v[2]; *target = v[3]; *disp = v[4];
  *fresh = v[5];
  return true;
}

/**
## Boundary conditions

Every velocity boundary is shifted by $-U$. `u.n[embed]` is the axial
component on the embedded wall and `u.t[embed]` the radial one (Basilisk
convention, see `src/test/couette.c`). */

u.n[left]  = dirichlet (POISEUILLE_FRAME);
u.t[left]  = dirichlet (0.);
p[left]    = neumann (0.);
pf[left]   = neumann (0.);
f[left]    = dirichlet (0.);

u.n[right] = dirichlet (POISEUILLE_FRAME);
u.t[right] = dirichlet (0.);
p[right]   = neumann (0.);
pf[right]  = neumann (0.);
f[right]   = dirichlet (0.);

u.n[embed] = dirichlet (-frame.U);
u.t[embed] = dirichlet (0.);

int main (int argc, char const *argv[])
{
  params_init_from_argv (argc, argv);

  CaseNo   = param_int ("CaseNo", 1000);
  MAXlevel = param_int ("MAXlevel", 10);
  MINlevel = param_int ("MINlevel", 4);
  padLevel = param_int ("padLevel", MAXlevel - 2);

  Ca      = param_double ("Ca", 0.05);
  CaPrev  = param_double ("CaPrev", 0.);
  La      = param_double ("La", 1.);
  muR     = param_double ("muR", 1e-2);
  rhoR    = param_double ("rhoR", 1e-3);
  Rtube   = param_double ("Rtube", 0.7);
  Rb0frac = param_double ("Rb0frac", 0.8);
  xRear   = param_double ("xRear", 1.05);
  Ldomain = param_double ("Ldomain", 8.);

  bTol     = param_double ("bTol", 2e-3);
  convHold = param_int    ("convHold", 3);
  wRenew   = param_double ("wRenew", 0.5);
  burnRenew = param_double ("burnRenew", 1.0);
  posTol   = param_double ("posTol", 0.02);

  filmBins      = param_int    ("filmBins", 128);
  filmCoverage  = param_double ("filmCoverage", .95);
  filmFlatTol   = param_double ("filmFlatTol", .05);
  freshFracMin  = param_double ("freshFracMin", .98);
  speedTol      = param_double ("speedTol", .02);
  shapeTol      = param_double ("shapeTol", .02);
  filmCells     = param_double ("filmCells", 4.);
  filmMinLevel  = param_int    ("filmMinLevel", 0);

  uRel   = param_double ("uRel", 1e-2);
  dRel   = param_double ("dRel", 1e-2);
  csErr  = param_double ("csErr", 1e-2);
  VelErr = uRel*Ca;
  DErr   = dRel*Ca/Rtube;

  travelR = param_double ("travelR", 8.);
  tmax    = param_double ("tmax", travelR/Ca);
  tsnap   = param_double ("tsnap", tmax/200.);
  tRamp   = param_double ("tRamp", 1.);
  filmSampleDt = param_double ("filmSampleDt", min(tsnap, 0.25/(8.*Ca)));

  DT = param_double ("dtmax", 1e-2);
  solverTol = param_double ("solverTol", 1e-4);

  tau = param_double ("tau", 1.);
  Uframe0 = param_double ("Uframe0", 0.);
  prescribedU = param_bool ("prescribedU", false);
  const char * outletKey = param_string ("outlet", "poiseuille");
  if (!strcmp (outletKey, "poiseuille")) outletMode = 0;
  else if (!strcmp (outletKey, "open-rear")) outletMode = 1;
  else if (!strcmp (outletKey, "parent")) outletMode = 2;
  else {
    fprintf (ferr, "ERROR: outlet must be poiseuille, open-rear or parent.\n");
    return 1;
  }

  Rb0  = Rb0frac*Rtube;
  Lcyl = 4.*(1. - cube(Rb0))/(3.*sq(Rb0));
  Xb0  = xRear + Rb0 + Lcyl/2.;
  vol0 = 4.*pi/3.;
  freshFront = param_double ("freshFront", Xb0 + Lcyl/2. + Rb0);
  xTarget = param_double ("xTarget", Xb0);

  const double finiteParameters[] = {
    Ca, CaPrev, La, muR, rhoR, Rtube, Rb0frac, xRear, Ldomain, travelR,
    tmax, tsnap, DT, tRamp, bTol, uRel, dRel, csErr, filmCoverage,
    filmFlatTol, freshFracMin, speedTol, shapeTol, filmCells, filmSampleDt,
    solverTol, freshFront, tau, Uframe0, xTarget, wRenew, burnRenew, posTol
  };
  for (unsigned int j = 0; j < sizeof(finiteParameters)/sizeof(double); j++)
    if (!isfinite(finiteParameters[j])) {
      fprintf (ferr, "ERROR: Runtime parameters must be finite.\n");
      return 1;
    }
  if (CaseNo < 1000 || MAXlevel <= 0 || MINlevel <= 0 ||
      MINlevel > MAXlevel || padLevel < MINlevel || padLevel > MAXlevel ||
      Ca <= 0. || CaPrev < 0. || La <= 0. || muR <= 0. ||
      rhoR <= 0. || Rtube <= 0. || Rtube >= 1. ||
      Rb0frac <= 0. || Rb0frac >= 1. || xRear <= 0. || Ldomain <= 0. ||
      tmax <= 0. || tsnap <= 0. || DT <= 0. || tRamp < 0. ||
      travelR <= 0. || bTol <= 0. || convHold < 1 || uRel <= 0. ||
      dRel <= 0. || csErr <= 0. || filmBins < 32 ||
      filmBins > CENTRAL_FILM_MAX_BINS || filmCoverage <= 0. ||
      filmCoverage > 1. || filmFlatTol <= 0. || freshFracMin < 0. ||
      freshFracMin > 1. || speedTol <= 0. || shapeTol <= 0. ||
      filmCells < 0. || filmMinLevel < 0 || filmMinLevel > MAXlevel ||
      filmSampleDt <= 0. || solverTol <= 0. || tau <= 0. || wRenew <= 0. ||
      burnRenew < 0. || posTol <= 0.) {
    fprintf (ferr, "ERROR: Invalid runtime parameters.\n");
    return 1;
  }
  if (Xb0 + Lcyl/2. + Rb0 >= Ldomain - 2.*Rtube) {
    fprintf (ferr, "ERROR: initial capsule leaves less than two tube radii of "
             "front pad; increase Ldomain or reduce xRear.\n");
    return 1;
  }

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
  sprintf (frameFile, "frame-state");

  int restartPresent = 0;
  if (pid() == 0) {
    struct stat restartStatus;
    restartPresent = stat(dumpFile, &restartStatus) == 0 &&
                     restartStatus.st_size > 0;
  }
  int logSetupError = 0;
  if (pid() == 0 && !prepare_case_log()) {
    fprintf (ferr, "HARDFAIL_LOG_SCHEMA: existing log %s is incompatible with "
             "this case or cannot be read.\n", logFile);
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

  rho1 = La*rhoR; mu1 = muR;
  rho2 = La;      mu2 = 1.;
  f.sigma = 1.;
  a = av;

  RenewalConfig rc = {
    .window_renewals = wRenew, .burn_renewals = burnRenew,
    .film_tol = bTol, .speed_tol = speedTol, .position_tol = posTol*Rtube,
    .length_tol = shapeTol, .minimum_coverage = filmCoverage,
    .maximum_spread = filmFlatTol, .minimum_fresh = freshFracMin,
    .minimum_cells = filmCells, .hold_windows = convHold
  };
  renewal_observer_init (&renewal, rc);
  comoving_frame_init (&frame, Uframe0, xTarget, tau, prescribedU);

  TOLERANCE = solverTol;
  CFL = 0.5;

  if (pid() == 0) {
    fprintf (ferr, "CaseNo=%d MAXlevel=%d MINlevel=%d padLevel=%d Ca=%g "
             "CaPrev=%g La=%g muR=%g rhoR=%g Rtube=%g Rb0=%g Lcyl=%g L0=%g "
             "tmax=%g solverTol=%g outlet=%s tau=%g Uframe0=%g xTarget=%g "
             "prescribedU=%d\n", CaseNo, MAXlevel, MINlevel, padLevel, Ca,
             CaPrev, La, muR, rhoR, Rtube, Rb0, Lcyl, Ldomain, tmax,
             solverTol, outletKey, tau, Uframe0, xTarget, prescribedU);
    fprintf (ferr, "Renewal windows: wRenew=%g burnRenew=%g convHold=%d "
             "bTol=%g speedTol=%g posTol=%g\n", wRenew, burnRenew, convHold,
             bTol, speedTol, posTol);
    fprintf (ferr, "Logging to %s\n", logFile);
  }
  run();
  return runFailed ? 1 : 0;
}

/** Open (traction-free) ends are selected at run time by overriding the
Dirichlet defaults declared above. */
event defaults (i = 0)
{
  if (outletMode == 1) {          // open rear end
    u.n[left] = neumann (0.);
    u.t[left] = neumann (0.);
    p[left]   = dirichlet (0.);
    pf[left]  = dirichlet (0.);
  }
  else if (outletMode == 2) {     // parent: open front end
    u.n[right] = neumann (0.);
    u.t[right] = neumann (0.);
    p[right]   = dirichlet (0.);
    pf[right]  = dirichlet (0.);
  }
}

event init (t = 0)
{
  restoredRun = restore (file = dumpFile);
  if (!restoredRun && restartFilePresent) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_RESTORE: restart exists but could not be "
               "restored.\n");
    terminate_init_failure();
  }
  if (!restoredRun && logHasContent) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_LOG_REUSE: %s exists but no restart was "
               "restored.\n", logFile);
    terminate_init_failure();
  }
  if (!restoredRun) {
    refine (y < 1.05*Rtube && level < MAXlevel - 2);
    refine (y < 1.05*Rtube &&
            x > Xb0 - Lcyl/2. - Rb0 - 2.*Rtube &&
            x < Xb0 + Lcyl/2. + Rb0 + 2.*Rtube && level < MAXlevel);
    tube_solid (Rtube);
    fraction (f, Rb0 - sqrt (sq (max (fabs (x - Xb0) - Lcyl/2., 0.))
                             + sq(y)));
    vof_solid_cleanup (f);
    tStart = 0.;
  }
  else {
    if (!close_parameter(L0, Ldomain)) {
      if (pid() == 0)
        fprintf (ferr, "HARDFAIL_DOMAIN_MISMATCH: restart L0=%g differs from "
                 "Ldomain=%g.\n", L0, Ldomain);
      terminate_init_failure();
    }
    tube_solid (Rtube);
#if TREE && AXI
    foreach_cell()
      if (cs[] >= 1.)
        cm[] = y;
#endif
    vof_solid_cleanup (f);
    // The global clock is not yet the restored time inside this event, so the
    // ramp origin is taken from the first time step (see frameControl).
    tStart = -1.;
    double tt, U, dUdt, target, disp, fresh;
    // A frame-state file belongs to this case directory: continuation
    // stations never copy it, so its presence identifies a same-case resume.
    if (read_frame_state (&tt, &U, &dUdt, &target, &disp, &fresh)) {
      // same-case resume: keep the frame exactly where the dump left it
      comoving_frame_init (&frame, U, target, tau, prescribedU);
      frame.dUdt = dUdt;
      frame.displacement = disp;
      freshFront = fresh;
      CaPrev = Ca;          // no new ramp on a same-case resume
      if (pid() == 0)
        fprintf (ferr, "# same-case resume: t=%g U=%g dUdt=%g target=%g "
                 "displacement=%g freshFront=%g\n", tt, U, dUdt, target, disp,
                 fresh);
    }
    else if (pid() == 0)
      fprintf (ferr, "# continuation seed: Uframe0=%g xTarget=%g CaPrev=%g "
               "freshFront=%g\n", Uframe0, xTarget, CaPrev, freshFront);
  }

  double front, rear, minimumFilm, centroid, axisTip;
  int components;
  if (!current_drop_geometry (&front, &rear, &minimumFilm, &centroid,
                              &axisTip, &components)) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_GEOMETRY: cannot reconstruct the dispersed "
               "phase after initialisation.\n");
    terminate_init_failure();
  }
  if (!restoredRun && param_string ("xTarget", NULL) == NULL) {
    frame.target = centroid;
    xTarget = centroid;
  }
  bool initialFilmOK = measure_central_film
    (rear, front, centroid, &latestFilmMeasurement);
  if (initialFilmOK)
    latestFilmValid = true;
  refine_measured_film();
  embed_axi_metric_sync();
  vof_solid_cleanup (f);
  if (!restoredRun) {
    // start from the frame-consistent fully developed profile in the liquid
    foreach()
      u.x[] = cs[] > 0. && f[] < 0.5 ?
        2.*CaPrev*(1. - sq(y/Rtube)) - frame.U : 0.;
  }
  if (pid() == 0)
    fprintf (ferr, "# init: front=%g rear=%g centroid=%g target=%g U=%g "
             "components=%d\n", front, rear, centroid, frame.target, frame.U,
             components);
}

/** Frame controller: runs before the solver's `last` events each step. */
event frameControl (i++)
{
  if (tStart < 0.)
    tStart = t;   // first step after a restore: the clock is now restored
  double xc = bubble_centroid();
  comoving_frame_update (&frame, xc, t);
}

event acceleration (i++)
{
  double ax = -frame.dUdt;
  foreach_face(x)
    av.x[] = ax;
  foreach_face(y)
    av.y[] = 0.;
}

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
  if (nonfinite || !isfinite(ke) || !isfinite(frame.U) ||
      !isfinite(frame.dUdt)) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_NONFINITE: nonfinite solver state at i=%d "
               "t=%g dt=%g U=%g dUdt=%g.\n", i, t, dt, frame.U, frame.dUdt);
    runFailed = 1;
    return 1;
  }
  if (ke > 1e3 && i > 10) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_ENERGY: kinetic energy %.6e exceeds 1e3 at "
               "i=%d t=%g.\n", ke, i, t);
    runFailed = 1;
    dump (file = "hardfail-energy");
    return 1;
  }
}

event adapt (i++)
{
  scalar KAPPA[];
  curvature (f, KAPPA);
  scalar Dmag[];
  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    double D22 = (y > 1e-10) ? u.y[]/y : D11;
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

event writingFiles (t = 0; t += tsnap; t <= tmax)
{
  dump (file = dumpFile);
  write_frame_state();
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

/** JSON has no NaN: non-finite values are written as null. */
static void json_number (FILE * fp, const char * key, double value,
                         const char * tail)
{
  if (isfinite (value))
    fprintf (fp, "  \"%s\": %.17g%s\n", key, value, tail);
  else
    fprintf (fp, "  \"%s\": null%s\n", key, tail);
}

/** Writes the station receipt and the interface polyline. */
static void write_station (const char * status, RenewalWindow * w)
{
  if (pid() != 0)
    return;
  FILE * fp = fopen ("station.json", "w");
  if (!fp)
    return;
  double h = w ? w->film : latestFilmMeasurement.film;
  double U = w ? w->U : frame.U;
  // Basilisk traps floating-point exceptions: guard the derived ratios.
  double stagnant = Rtube > h ? sq(Rtube/(Rtube - h)) : NAN;
  double speedRatio = Ca > 0. ? U/Ca : NAN;
  double stagnantResidual = (U > 0. && stagnant > 0.) ?
    speedRatio/stagnant - 1. : NAN;
  double mobility = U > 0. ? 1. - Ca/U : NAN;
  fprintf (fp, "{\n  \"status\": \"%s\",\n  \"CaseNo\": %d,\n"
           "  \"MAXlevel\": %d,\n  \"padLevel\": %d,\n  \"outlet\": %d,\n"
           "  \"i\": %d,\n  \"hold_count\": %d,\n", status, CaseNo, MAXlevel,
           padLevel, outletMode, iter, w ? w->hold_count : 0);
  json_number (fp, "Ca_in", Ca, ",");
  json_number (fp, "La", La, ",");
  json_number (fp, "muR", muR, ",");
  json_number (fp, "rhoR", rhoR, ",");
  json_number (fp, "Rtube", Rtube, ",");
  json_number (fp, "Ldomain", Ldomain, ",");
  json_number (fp, "tau", tau, ",");
  json_number (fp, "t", t, ",");
  json_number (fp, "U", frame.U, ",");
  json_number (fp, "Ca_b", U, ",");
  json_number (fp, "v_rel", frame.v_rel, ",");
  json_number (fp, "dUdt", frame.dUdt, ",");
  json_number (fp, "xTarget", frame.target, ",");
  json_number (fp, "displacement", frame.displacement, ",");
  json_number (fp, "freshFront", freshFront, ",");
  json_number (fp, "h", h, ",");
  json_number (fp, "h_over_Rtube", h/Rtube, ",");
  json_number (fp, "h_over_R0", h, ",");
  json_number (fp, "speed_ratio", speedRatio, ",");
  json_number (fp, "stagnant_film_residual", stagnantResidual, ",");
  json_number (fp, "m", mobility, "");
  fprintf (fp, "}\n");
  fclose (fp);
  FILE * fi = fopen ("interface.dat", "w");
  if (fi) {
    output_facets (f, fi);
    fclose (fi);
  }
}

event logWriting (i++)
{
  if (t + 1e-12 < nextFilmSampleTime)
    return 0;
  nextFilmSampleTime = t + filmSampleDt;

  double ke = latestSafetyKE;
  double vol = 2.*pi*statsf(f).sum;
  double xTipF, xTipR, bFilm, centroid, axisTip;
  int components;
  bool geometryOK = current_drop_geometry
    (&xTipF, &xTipR, &bFilm, &centroid, &axisTip, &components);
  if (!geometryOK) {
    if (pid() == 0)
      fprintf (ferr, "HARDFAIL_GEOMETRY: no dispersed-phase component at "
               "t=%g.\n", t);
    runFailed = 1;
    return 1;
  }
  // film deposited since the last Ca change has moved rearward by the
  // frame displacement
  freshFront = xTipF - frame.displacement;
  bool centralOK = measure_central_film
    (xTipR, xTipF, centroid, &latestFilmMeasurement);
  if (centralOK)
    latestFilmValid = true;
  double central = centralOK ? latestFilmMeasurement.film : NAN;
  double spread = centralOK ? latestFilmMeasurement.spatial_spread : NAN;
  double coverage = centralOK ? latestFilmMeasurement.coverage : 0.;
  double fresh = centralOK ? latestFilmMeasurement.fresh_fraction : 0.;
  double cells = centralOK ? latestFilmMeasurement.minimum_cells : 0.;
  double Qrear, Qfront, pFront;
  boundary_fluxes (&Qrear, &Qfront, &pFront);
  double pB = bubble_pressure();
  double pExcess = pB - pFront - 8.*ca_now()*(L0 - xTipF)/sq(Rtube);
  double delta_tail = axisTip - xTipR;

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
        fprintf (fp, "# i dt t ke dVol/Vol0 xTipF xTipR bFilm bCentral "
                 "spatialSpread coverage freshFraction minCells centroid "
                 "length Uframe dUdt vRel xcErr Qrear Qfront pBubble "
                 "pFront pExcess deltaTail components CaNow\n");
        logNeedsHeader = false;
      }
      fprintf (fp, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e "
               "%.6e %.6e %.6e %.6e %.9e %.6e %.6e %.6e %.9e %.9e %.6e "
               "%.6e %.6e %.6e %d %.6e\n",
               i, dt, t, ke, (vol - vol0)/vol0, xTipF, xTipR, bFilm,
               central, spread, coverage, fresh, cells, centroid,
               xTipF - xTipR, frame.U, frame.dUdt, frame.v_rel, frame.error,
               Qrear, Qfront, pB, pFront, pExcess, delta_tail, components,
               ca_now());
      fclose (fp);
    }
    fprintf (ferr, "%d %.4e %.5e U=%.6e dUdt=%.2e vRel=%.2e err=%.2e "
             "h=%.5e/%.5e spread=%.2e fresh=%.2f cells=%.1f Q=%.3e/%.3e "
             "dV=%.2e comp=%d\n", i, dt, t, frame.U, frame.dUdt, frame.v_rel,
             frame.error, central, central/Rtube, spread, fresh, cells,
             Qrear, Qfront, (vol - vol0)/vol0, components);
  }
  if (logError) {
#if _MPI
    MPI_Abort (MPI_COMM_WORLD, 1);
#endif
    runFailed = 1;
    return 1;
  }

  if (bFilm < 4.*L0/(1 << MAXlevel) && i > 10 && pid() == 0)
    fprintf (ferr, "WARNING: film resolved by fewer than 4 cells at t=%g.\n", t);
  if (xTipF > L0 - 1.5*Rtube || xTipR < 0.5*Rtube) {
    if (pid() == 0)
      fprintf (ferr, "INCOMPLETE_DRIFT: bubble left the window at t=%g "
               "(front=%g rear=%g).\n", t, xTipF, xTipR);
    runFailed = 1;
    dump (file = dumpFile);
    write_frame_state();
    write_station ("INCOMPLETE_DRIFT", NULL);
    return 1;
  }

  if (centralOK && tStart >= 0. && t - tStart >= tRamp && frame.U > 0.) {
    RenewalSample s = {
      .time = t, .film = latestFilmMeasurement.film, .U = frame.U,
      .xc = centroid, .length = xTipF - xTipR, .error = frame.error,
      .spread = spread, .coverage = coverage, .fresh = fresh, .cells = cells
    };
    RenewalWindow w;
    double renewalTime = (xTipF - xTipR)/frame.U;
    if (renewal_observer_push (&renewal, s, renewalTime, &w)) {
      if (pid() == 0)
        fprintf (ferr, "# renewal window: t=[%.3f,%.3f] h=%.6e (h/Rt=%.6e) "
                 "U=%.6e filmDrift=%.3e speedDrift=%.3e lengthDrift=%.3e "
                 "errMax=%.3e spread=%.3e fresh=%.3f minCells=%.2f "
                 "hold=%d/%d\n", w.t_start, w.t_end, w.film, w.film/Rtube,
                 w.U, w.film_drift, w.speed_drift, w.length_drift,
                 w.error_abs_max, w.max_spread, w.min_fresh, w.min_cells,
                 w.hold_count, convHold);
      if (w.converged) {
        if (pid() == 0)
          fprintf (ferr, "SUCCESS: renewal-window stationarity at t=%g: "
                   "h=%.6e h/Rtube=%.6e U=Ca_b=%.6e Ca_in=%g\n", t, w.film,
                   w.film/Rtube, w.U, Ca);
        dump (file = dumpFile);
        write_frame_state();
        write_station ("SUCCESS", &w);
        return 1;
      }
    }
  }
}

event stopSimulation (t = tmax)
{
  runFailed = 1;
  if (pid() == 0)
    fprintf (ferr, "INCOMPLETE_TMAX: case %d reached tmax without "
             "stationarity: Ca %g.\n", CaseNo, Ca);
  dump (file = dumpFile);
  write_frame_state();
  write_station ("INCOMPLETE_TMAX", NULL);
  return 1;
}
