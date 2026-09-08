/**
# Axisymmetric embedded-tube geometry verification

Checks the exact volume and axial area of a straight tube on a deliberately
nonuniform, repeatedly adapted quadtree. No flow equations are solved.
The geometric integrals omit the common azimuthal factor of $2\pi$.
*/
#include "grid/quadtree.h"
#include "embed.h"
#include "axi.h"
#include "run.h"
#include "embed-vof-tube.h"

static const double RT = 0.7;
static int MAXLEVEL = 9, failed = 0;

int main (int argc, char ** argv)
{
  if (argc > 1)
    MAXLEVEL = atoi (argv[1]);
  if (MAXLEVEL < 6 || MAXLEVEL > 13) {
    fprintf (stderr, "FAIL invalid maximum level %d\n", MAXLEVEL);
    return 2;
  }
  size (16.);
  origin (0., 0.);
  init_grid (1 << 4);
  run();
  return failed;
}

event init (i = 0)
{
  /* Cross the wall with several resolution transitions from the outset. */
  refine (y < 1.1 && level < 6);
  refine (y < 0.9 && x > 0.75 && x < 3.25 && level < MAXLEVEL);
  refine (fabs(y - RT) < 0.08 && x > 6. && x < 10. &&
          level < MAXLEVEL - 1);
  tube_solid (RT);
}

event geometry_cycle (i = 0; i <= 8; i++)
{
  /* Moving compact features force both refinement and coarsening through the
     wall and through changing MPI ownership/halo topology. */
  scalar marker[];
  double xc = 1.5 + 1.65*i;
  foreach()
    marker[] = exp (-sq((x - xc)/0.35) - sq((y - RT)/0.055));
  adapt_wavelet ({marker, cs}, (double[]){2e-2, 1e-3},
                 MAXLEVEL, 4);
  embed_axi_metric_sync();

  /* Synchronise fractions and metrics before evaluating the geometry. */
  boundary ({cs, fs, cm, fm});

  double volume = 0., inlet = 0., outlet = 0., divmax = 0.;
  int bad = 0;
  foreach (reduction(+:volume) reduction(max:divmax)) {
    volume += cm[]*sq(Delta);
    if (cs[] > 0.) {
      double residual = fabs (fm.x[1] - fm.x[]);
      if (residual > divmax)
        divmax = residual;
    }
  }
  foreach_boundary (left, reduction(+:inlet))
    inlet += fm.x[]*Delta;
  foreach_boundary (right, reduction(+:outlet))
    outlet += fm.x[1]*Delta;
  foreach_face (reduction(+:bad))
    if ((is_active(cell) || is_active(neighbor(-1))) &&
        fs.x[] > 0. && fm.x[] <= 0.)
      bad++;

  const double exact_volume = L0*sq(RT)/2.;
  const double exact_area = sq(RT)/2.;
  double ev = fabs(volume - exact_volume);
  double ei = fabs(inlet - exact_area);
  double eo = fabs(outlet - exact_area);
  int pass = ev <= 5e-12 && ei <= 5e-12 && eo <= 5e-12 &&
             divmax <= 5e-12 && bad == 0;

  if (pid() == 0) {
    if (i == 0)
      fprintf (stdout,
               "cycle,volume,volume_error,inlet_area,inlet_error,"
               "outlet_area,outlet_error,constant_flux_divmax,bad_faces,status\n");
    fprintf (stdout, "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
             "%.17g,%d,%s\n", i, volume, ev, inlet, ei, outlet, eo,
             divmax, bad, pass ? "PASS" : "FAIL");
    fflush (stdout);
  }
  if (!pass) {
    if (pid() == 0)
      fprintf (stderr, "FAIL geometry predicate at cycle %d\n", i);
    failed = 1;
    return 1;
  }
  if (i == 8) {
    if (pid() == 0)
      fprintf (stdout, "PASS\n");
    return 1;
  }
}
