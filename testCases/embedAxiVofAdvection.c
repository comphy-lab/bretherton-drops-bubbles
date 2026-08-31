/**
# embedAxiVofAdvection.c

**Verification** of the embed + axi + VOF coupling against an exact
solution (evidence source: analytical solution of the advection
equation, so this is verification, not a software smoke test).

Adapted from the upstream compatibility test
[`src/test/missing_metric.c`](http://basilisk.fr/src/test/missing_metric.c),
which guards the historical incompatibility between embedded boundaries,
the axisymmetric metric and VOF advection: with both `AXI` and `EMBED`,
the metric is $f_m = y f_s$, and VOF fluxes must include it. Two changes
are made here:

1. a *genuine* solid wall is embedded at $y = y_w$ (the upstream test
   sets `cs = 1` everywhere, exercising the code path but no real
   geometry);
2. the test is self-verifying: it compares the interface position with
   the exact solution and prints `PASS`/`FAIL`.

An axisymmetric interface initially at radius $r_0$ is advected by the
divergence-free radial field $u_r = 1/r$, so its exact position is
$r(t) = \sqrt{r_0^2 + 2t}$. The wall at $y_w = 1.9$ is not reached
within the simulated time ($r(0.8) \simeq 1.74$), so the exact solution
is unaffected; the embedded boundary must simply not corrupt the
advection.
*/

#include "embed.h"
#include "axi.h"
#define ro 1.2
#define veloc(r) (1./(r))
#define fpos(t) sqrt(sq(ro) + 2.*(t))
#define ywall 1.9
#define TOL_REL 1e-2

#include "advection.h"
#include "vof.h"
#include "curvature.h"
#include "embed-vof-tube.h"

scalar f[], * tracers = NULL, * interfaces = {f};

double maxerr = 0.;

int main()
{
  Y0 = 1.0;
  L0 = 1.0;
  DT = HUGE;
  N = 64;
  run();
}

event init (i = 0)
{
  solid (cs, fs, ywall - y);
  embed_axi_metric_sync();

  fraction (f, ro - y);
  vof_solid_cleanup (f);

  /**
  The face velocity carries the full metric $f_m = y f_s$, so fluxes
  through embedded faces vanish identically. */

  foreach_face(y)
    u.y[] = veloc(y)*fm.y[];
}

/**
## Results

The maximum interface radius is compared with the exact solution at
regular intervals.
*/
event prof_pos (i += 5; t <= 0.8)
{
  scalar pos[];
  position (f, pos, {0, 1.});
  double rnum = statsf(pos).max, rex = fpos(t);
  double err = fabs (rnum - rex)/rex;
  if (err > maxerr)
    maxerr = err;
  fprintf (stderr, "%g %g %g %g\n", t, rnum, rex, err);
}

event end (t = 0.8)
{
  printf ("max relative interface-position error: %g\n", maxerr);
  printf (maxerr < TOL_REL ? "PASS\n" : "FAIL\n");
}
