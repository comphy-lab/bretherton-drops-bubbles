/**
# Restart regridding contract for the Bretherton tube solver

This is a bounded software test of the transfer operators used when a saved
production state is refined by one or two levels.  It restores an actual
Basilisk dump, rebuilds the same embedded axisymmetric tube, applies the same
film-floor geometry as `bretherton.c`, and stops before time integration.

The test reports axisymmetric dispersed volume and centroid, mixture momentum,
kinetic energy, exact straight-tube geometry errors, solid contamination and
the discrete leaf face-flux divergence before and after one projection.  On an
adaptive embedded tree this direct leaf readback is a diagnostic, not the
multilevel norm controlled by the Poisson solver.  The projection uses the
operator-normalised value `dt = 1`; it is not a substitute for a normal
production timestep with the complete solver event sequence and does not by
itself establish pointwise incompressibility.
Running the same command with serial and MPI builds supplies the parallel
agreement check.  These measurements establish only the restart/regrid
implementation contract; they do not establish grid convergence or physical
stationarity.

~~~bash
qcc -I../src-local -O2 -Wall -disable-dimensions regridTube.c -o regridTube -lm
./regridTube SNAPSHOT LEVELS RTUBE FILM CA RHO1 RHO2 [TOLERANCE]
~~~

`LEVELS` must be one or two. `FILM` is the dimensional central-film thickness
used by the source checkpoint (for example `Rtube*(b/Rtube)`).  The remaining
arguments reproduce the production boundary conditions and density law.
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "embed-vof-tube.h"
#include <errno.h>
#include <float.h>

static const char * snapshotFile;
static int levelIncrement, sourceLevel, targetLevel;
static double Rtube, filmThickness, Ca, tRamp = 1.;
static int failed = 0;

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

typedef struct {
  long leaves, mixed, cut;
  int level_min, level_max, mixed_level_min, mixed_level_max;
  double volume, centroid_x, centroid_y;
  double momentum_x, momentum_y, kinetic_energy;
  double tube_volume_error, inlet_area_error, outlet_area_error;
  double cs_error, fs_error, cm_error, fm_error, solid_fraction_max;
  double leaf_flux_divergence_diagnostic;
  double physical_divergence_full_max, physical_divergence_cut_max;
  double minimum_positive_cm;
} RegridMetrics;

static bool parse_double (const char * text, double * value)
{
  errno = 0;
  char * end = NULL;
  *value = strtod (text, &end);
  return errno == 0 && end != text && *end == '\0' && isfinite (*value);
}

static bool parse_increment (const char * text, int * value)
{
  errno = 0;
  char * end = NULL;
  long parsed = strtol (text, &end, 10);
  if (errno || end == text || *end != '\0' || (parsed != 1 && parsed != 2))
    return false;
  *value = (int) parsed;
  return true;
}

/** Reconstruct the face flux exactly as `centered.h` does after user init. */
static void rebuild_face_flux (void)
{
  event ("properties");
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  boundary ((scalar *){uf});
}

static void divergence_metrics (RegridMetrics * metrics)
{
  double rawMaximum = 0., fullMaximum = 0., cutMaximum = 0., minimumCm = HUGE;
  foreach (serial, noauto)
    if (cm[] > 0.) {
      double divergence = 0.;
      foreach_dimension()
        divergence += uf.x[1] - uf.x[];
      const double raw = fabs (divergence/Delta);
      const double physical = raw/cm[];
      rawMaximum = max (rawMaximum, raw);
      minimumCm = min (minimumCm, cm[]);
      if (cs[] >= 1.)
        fullMaximum = max (fullMaximum, physical);
      else if (cs[] > 0.)
        cutMaximum = max (cutMaximum, physical);
    }
#if _MPI
  mpi_all_reduce (rawMaximum, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (fullMaximum, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (cutMaximum, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (minimumCm, MPI_DOUBLE, MPI_MIN);
#endif
  metrics->leaf_flux_divergence_diagnostic = rawMaximum;
  metrics->physical_divergence_full_max = fullMaximum;
  metrics->physical_divergence_cut_max = cutMaximum;
  metrics->minimum_positive_cm = minimumCm < HUGE ? minimumCm : 0.;
}

/** Compare the live tube fractions and metrics with the exact straight tube. */
static void geometry_errors (RegridMetrics * metrics)
{
  double tubeVolume = 0., inlet = 0., outlet = 0.;
  double ecs = 0., efs = 0., ecm = 0., efm = 0., fsolid = 0.;

  scalar cmv = cm;
  face vector fmv = fm;
  foreach (reduction(+:tubeVolume) reduction(max:ecs)
           reduction(max:ecm) reduction(max:fsolid)) {
    tubeVolume += dv();
    double csExact = clamp ((Rtube - y)/Delta + .5, 0., 1.);
    double cmExact = tube_axial_face_metric (y, Delta);
    ecs = max (ecs, fabs (cs[] - csExact));
    ecm = max (ecm, fabs (cmv[] - cmExact));
    if (cs[] <= 0.)
      fsolid = max (fsolid, fabs (f[]));
  }
  foreach_face (x, reduction(max:efs) reduction(max:efm)) {
    double exact = clamp ((Rtube - y)/Delta + .5, 0., 1.);
    efs = max (efs, fabs (fs.x[] - exact));
    efm = max (efm, fabs (fmv.x[] - tube_axial_face_metric(y, Delta)));
  }
  foreach_face (y, reduction(max:efs) reduction(max:efm)) {
    double fsExact = y < Rtube ? 1. : 0.;
    double fmExact = fsExact ? max(y, 1e-20) : 0.;
    efs = max (efs, fabs (fs.y[] - fsExact));
    efm = max (efm, fabs (fmv.y[] - fmExact));
  }
  foreach_boundary (left, reduction(+:inlet))
    inlet += fm.x[]*Delta;
  foreach_boundary (right, reduction(+:outlet))
    outlet += fm.x[1]*Delta;

  const double exactVolume = L0*sq(Rtube)/2.;
  const double exactArea = sq(Rtube)/2.;
  metrics->tube_volume_error = fabs (tubeVolume - exactVolume);
  metrics->inlet_area_error = fabs (inlet - exactArea);
  metrics->outlet_area_error = fabs (outlet - exactArea);
  metrics->cs_error = ecs;
  metrics->fs_error = efs;
  metrics->cm_error = ecm;
  metrics->fm_error = efm;
  metrics->solid_fraction_max = fsolid;
}

static RegridMetrics measure_state (void)
{
  RegridMetrics metrics = {.level_min = 30, .mixed_level_min = 30};
  double volume = 0., mx = 0., my = 0., px = 0., py = 0., ke = 0.;
  long leaves = 0, mixed = 0, cut = 0;
  int levelMin = 30, levelMax = 0, mixedMin = 30, mixedMax = 0;

  foreach (reduction(+:volume) reduction(+:mx) reduction(+:my)
           reduction(+:px) reduction(+:py) reduction(+:ke)
           reduction(+:leaves) reduction(+:mixed) reduction(+:cut)
           reduction(min:levelMin) reduction(max:levelMax)
           reduction(min:mixedMin) reduction(max:mixedMax)) {
    const double element = 2.*pi*dv();
    const double dispersed = f[]*element;
    const double density = rho(f[]);
    leaves++;
    levelMin = min (levelMin, level);
    levelMax = max (levelMax, level);
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      mixed++;
      mixedMin = min (mixedMin, level);
      mixedMax = max (mixedMax, level);
    }
    if (cs[] > 0. && cs[] < 1.)
      cut++;
    volume += dispersed;
    mx += x*dispersed;
    my += y*dispersed;
    px += density*u.x[]*element;
    py += density*u.y[]*element;
    ke += 0.5*density*(sq(u.x[]) + sq(u.y[]))*element;
  }

  metrics.leaves = leaves;
  metrics.mixed = mixed;
  metrics.cut = cut;
  metrics.level_min = levelMin;
  metrics.level_max = levelMax;
  metrics.mixed_level_min = mixed ? mixedMin : -1;
  metrics.mixed_level_max = mixed ? mixedMax : -1;
  metrics.volume = volume;
  metrics.centroid_x = volume > 0. ? mx/volume : nodata;
  metrics.centroid_y = volume > 0. ? my/volume : nodata;
  metrics.momentum_x = px;
  metrics.momentum_y = py;
  metrics.kinetic_energy = ke;
  divergence_metrics (&metrics);
  geometry_errors (&metrics);
  return metrics;
}

static void print_metrics (const char * stage, RegridMetrics m)
{
  if (pid() != 0)
    return;
  // iter is the global event counter declared by Basilisk's grid/events.h.
  printf ("%s,%d,%.17g,%d,%ld,%d,%d,%ld,%d,%d,%ld,"
          "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
          "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
          "%.17g,%.17g,%.17g\n",
          stage, npe(), t, iter, m.leaves, m.level_min, m.level_max,
          m.mixed, m.mixed_level_min, m.mixed_level_max, m.cut,
          m.volume, m.centroid_x, m.centroid_y,
          m.momentum_x, m.momentum_y, m.kinetic_energy,
          m.tube_volume_error, m.inlet_area_error, m.outlet_area_error,
          m.cs_error, m.fs_error, m.cm_error, m.fm_error,
          m.solid_fraction_max, m.leaf_flux_divergence_diagnostic,
          m.physical_divergence_full_max, m.physical_divergence_cut_max,
          m.minimum_positive_cm);
  fflush (stdout);
}

int main (int argc, char ** argv)
{
  if (argc != 8 && argc != 9) {
    fprintf (stderr, "usage: %s SNAPSHOT LEVELS RTUBE FILM CA RHO1 RHO2 "
             "[TOLERANCE]\n", argv[0]);
    return 2;
  }
  snapshotFile = argv[1];
  TOLERANCE = 1e-4;
  if (!parse_increment(argv[2], &levelIncrement) ||
      !parse_double(argv[3], &Rtube) || !parse_double(argv[4], &filmThickness) ||
      !parse_double(argv[5], &Ca) || !parse_double(argv[6], &rho1) ||
      !parse_double(argv[7], &rho2) ||
      (argc == 9 && !parse_double(argv[8], &TOLERANCE)) ||
      Rtube <= 0. || filmThickness <= 0. || Ca <= 0. ||
      rho1 <= 0. || rho2 <= 0. || TOLERANCE <= 0.) {
    fprintf (stderr, "FAIL invalid argument\n");
    return 2;
  }

  mu1 = mu2 = 1.;
  f.sigma = 1.;
  DT = 1e-2;
  init_grid (1 << 4);
  run();
  return failed;
}

event init (t = 0)
{
  if (!restore (file = snapshotFile)) {
    if (pid() == 0)
      fprintf (stderr, "FAIL cannot restore %s\n", snapshotFile);
    failed = 1;
    return 1;
  }

  int restoredMaximumLevel = 0;
  foreach (reduction(max:restoredMaximumLevel))
    restoredMaximumLevel = max (restoredMaximumLevel, level);
  sourceLevel = restoredMaximumLevel;
  targetLevel = sourceLevel + levelIncrement;
  if (targetLevel > 30) {
    if (pid() == 0)
      fprintf (stderr, "FAIL target level %d is unsupported\n", targetLevel);
    failed = 1;
    return 1;
  }

  tube_solid (Rtube);
  vof_solid_cleanup (f);
  boundary ({f, u, cs, fs, cm, fm});
  rebuild_face_flux();

  if (pid() == 0) {
    printf ("stage,ranks,t,iter,leaves,level_min,level_max,mixed,"
            "mixed_level_min,mixed_level_max,cut,volume,centroid_x,centroid_y,"
            "momentum_x,momentum_y,kinetic_energy,tube_volume_error,"
            "inlet_area_error,outlet_area_error,cs_error,fs_error,cm_error,"
            "fm_error,solid_fraction_max,leaf_flux_divergence_diagnostic,"
            "physical_divergence_full_max,physical_divergence_cut_max,"
            "minimum_positive_cm\n");
    fprintf (stderr, "# source_level=%d target_level=%d L0=%.17g "
             "Rtube=%.17g film=%.17g Ca=%.17g rho1=%.17g rho2=%.17g\n",
             sourceLevel, targetLevel, L0, Rtube, filmThickness, Ca,
             rho1, rho2);
  }
  print_metrics ("restored", measure_state());

  double front = -HUGE, rear = HUGE;
  scalar xpos[];
  position (f, xpos, {1, 0});
  front = statsf(xpos).max;
  rear = statsf(xpos).min;
  if (!(front > rear)) {
    if (pid() == 0)
      fprintf (stderr, "FAIL cannot reconstruct interface tips\n");
    failed = 1;
    return 1;
  }

  const double wall0 = rear - Rtube, wall1 = front + Rtube;
  const double band = 2.*filmThickness;
  tube_refinement_geometry_begin();
  refine ((level < targetLevel && cs[] > 0. &&
           f[] > 1e-6 && f[] < 1. - 1e-6) ||
          (level < targetLevel && cs[] > 0. &&
           x > wall0 && x < wall1 && y + Delta/2. > Rtube - band));
  tube_refinement_geometry_end();
  embed_axi_metric_sync();
  vof_solid_cleanup (f);
  boundary ({f, u, cs, fs, cm, fm});
  rebuild_face_flux();
  print_metrics ("regridded_pre_projection", measure_state());

  /**
  `dt = 1` makes the multilevel projection residual directly comparable with
  `TOLERANCE`. The direct leaf divergence and its physical value divided by
  `cm` remain useful AMR diagnostics, but are not the convergence norm used by
  Basilisk's tree Poisson operator.
  */
  const double projectionDt = 1.;
  mgstats projection = project (uf, p, alpha, projectionDt, 0);
  RegridMetrics projected = measure_state();
  print_metrics ("regridded_post_projection", projected);

  int achieved = 0;
  foreach (reduction(max:achieved))
    if (f[] > 1e-6 && f[] < 1. - 1e-6)
      achieved = max (achieved, level);
  if (pid() == 0)
    fprintf (stderr, "# standalone_projection dt=%g resb=%.17g resa=%.17g "
             "iterations=%d nrelax=%d target=%d mixed_max=%d "
             "post_leaf_flux_divergence=%.17g\n",
             projectionDt, projection.resb, projection.resa, projection.i,
             projection.nrelax, targetLevel, achieved,
             projected.leaf_flux_divergence_diagnostic);
  const double roundoffAllowance =
    64.*DBL_EPSILON*max (1., projection.resa);
  const bool fieldsFinite =
    isfinite (projected.volume) && isfinite (projected.centroid_x) &&
    isfinite (projected.centroid_y) && isfinite (projected.momentum_x) &&
    isfinite (projected.momentum_y) && isfinite (projected.kinetic_energy) &&
    isfinite (projected.leaf_flux_divergence_diagnostic);
  const bool geometryOK =
    projected.tube_volume_error <= 5e-12 &&
    projected.inlet_area_error <= 5e-12 &&
    projected.outlet_area_error <= 5e-12 &&
    projected.cs_error <= 5e-12 && projected.fs_error <= 5e-12 &&
    projected.cm_error <= 5e-12 && projected.fm_error <= 5e-12 &&
    projected.solid_fraction_max == 0.;
  const bool projectionOK =
    isfinite (projection.resb) && isfinite (projection.resa) &&
    projection.resa <= TOLERANCE + roundoffAllowance;
  if (achieved != targetLevel || !fieldsFinite || !geometryOK || !projectionOK) {
    if (pid() == 0)
      fprintf (stderr, "FAIL regrid/projection predicate: tolerance=%.17g "
               "roundoff_allowance=%.17g\n", TOLERANCE,
               roundoffAllowance);
    failed = 1;
  }
  return 1;
}
