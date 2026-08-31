/**
# laplaceEmbedTube.c

**Verification** of surface tension inside an embedded tube (evidence
source: exact Young–Laplace solution).

A static spherical droplet of radius $R_d$ sits on the axis of an
embedded tube of radius $R_{tube}$, surrounded by the continuous wetting
phase. The exact solution is zero velocity and a pressure jump
$\Delta p = 2\sigma/R_d$. Any velocity is a spurious (parasitic)
current; any deviation of the measured jump from $2\sigma/R_d$ measures
the curvature error of the height-function method operating near an
embedded boundary it knows nothing about.

This exercises the full production stack of `simulationCases/bretherton.c`
— embed + axi + centered + two-phase + tension — in a configuration with
a known solution. The droplet surface stays $R_{tube} - R_d$ away from
the wall, mimicking the wetting-film separation of the Bretherton
configuration.

Pass criteria (capillary units, $\sigma = \mu = 1$):

- relative pressure-jump error below 2%,
- maximum spurious velocity below $10^{-3}\,\sigma/\mu$ at $t = 1$.
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
#include "tension.h"
#include "embed-vof-tube.h"

#define Rtube 0.7
#define Rd 0.4
#define Xd 1.0
#define tend 1.0
#define TOL_DP 2e-2
#define TOL_U 1e-3

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

int main()
{
  L0 = 2.0;
  N = 128;

  rho1 = 1.; mu1 = 1.;
  rho2 = 1.; mu2 = 1.;
  f.sigma = 1.;

  TOLERANCE = 1e-5;
  DT = 1e-2;

  run();
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
checked at $t = t_{end}$. The pressure jump is measured between cell
averages well inside ($r < R_d/2$) and well outside ($r > 3R_d/2$, fluid
only) the droplet.
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
  fprintf (stderr, "%g %.6e %.6e\n", t, dp, umax);
}

event end (t = tend)
{
  double dp, umax;
  dp_umax (&dp, &umax);
  double dpex = 2./Rd;
  double edp = fabs (dp - dpex)/dpex;
  printf ("pressure jump: %.6e (exact %.6e, rel. error %.3e)\n",
          dp, dpex, edp);
  printf ("max spurious velocity: %.6e\n", umax);
  printf (edp < TOL_DP && umax < TOL_U ? "PASS\n" : "FAIL\n");
}
