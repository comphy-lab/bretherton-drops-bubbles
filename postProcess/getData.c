/**
# getData.c

Sample velocity magnitude and viscous dissipation from a `bretherton.c`
snapshot onto a uniform grid, for rendering.

The dispersed phase is `f = 1` and the continuous wetting phase is
`f = 0` (see `simulationCases/bretherton.c`), so the local viscosity is
$\mu(f) = f\mu_R + (1-f)$ in units of $\mu_c$. This is the opposite
assignment to the bursting-bubble post-processors, where the liquid
carries `f = 1`; taking their mapping unchanged would paint the bubble
interior with the continuous-phase viscosity and overstate dissipation
inside it by two orders of magnitude.

The tube wall is an embedded boundary, so cells inside the solid carry
no meaningful velocity. `cs` is emitted alongside the fields and the
renderer masks on it; nothing here guesses a value for the solid.

## Usage

```
./getData <dump> <xmin> <ymin> <xmax> <ymax> <ny> <muR>
```

Columns on stderr: `x y cs f log10(mu D:D) |u| u_x`.

`u_x` is the axial component. The mean of `|u|` inside a bubble
exceeds its translation speed whenever the interior recirculates, so a
plug-velocity check must use `u_x`, not the magnitude.

Dissipation is reported as $\log_{10}(\mu\,\mathbf{D}\!:\!\mathbf{D})$,
matching the sibling convention; the true viscous dissipation rate is
$2\mu\,\mathbf{D}\!:\!\mathbf{D}$, so this is a constant $\log_{10}2$
below it and the colourbar is labelled accordingly.
*/

#include "utils.h"
#include "output.h"

scalar f[], cs[];
vector u[];
scalar D2c[], vel[], ux[];

int main (int argc, char const *argv[])
{
  if (argc < 8) {
    fprintf (ferr, "usage: %s <dump> <xmin> <ymin> <xmax> <ymax> <ny> <muR>\n",
             argv[0]);
    return 1;
  }
  char filename[4096];
  snprintf (filename, sizeof(filename), "%s", argv[1]);
  double xmin = atof(argv[2]), ymin = atof(argv[3]);
  double xmax = atof(argv[4]), ymax = atof(argv[5]);
  int ny = atoi(argv[6]);
  double muR = atof(argv[7]);

  if (ny <= 0 || xmax <= xmin || ymax <= ymin) {
    fprintf (ferr, "ERROR: need xmax>xmin, ymax>ymin and ny>0\n");
    return 1;
  }

  restore (file = filename);

  /**
  Axisymmetric strain-rate invariant. `D22` is the hoop term
  $u_y/y$, guarded on the axis. */

  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    double D22 = (y > 1e-10) ? u.y[]/y : 0.;
    double D33 = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    double D13 = 0.5*((u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2.*Delta));
    double D2 = sq(D11) + sq(D22) + sq(D33) + 2.*sq(D13);
    double muLocal = clamp(f[],0.,1.)*muR + (1. - clamp(f[],0.,1.));
    double diss = muLocal*D2;
    D2c[] = diss > 0. ? log(diss)/log(10.) : -10.;
    vel[] = sqrt (sq(u.x[]) + sq(u.y[]));
    ux[] = u.x[];
  }

  double Deltay = (ymax - ymin)/(ny - 1);
  int nx = (int)((xmax - xmin)/Deltay) + 1;
  if (nx <= 0) {
    fprintf (ferr, "ERROR: computed nx <= 0\n");
    return 1;
  }
  double Deltax = (xmax - xmin)/(nx - 1);

  for (int i = 0; i < nx; i++) {
    double xp = xmin + i*Deltax;
    for (int j = 0; j < ny; j++) {
      double yp = ymin + j*Deltay;
      /**
      `interpolate()` returns `nodata` outside the domain and inside the
      embedded solid. Emitting that sentinel as a number would let a
      renderer plot 1e30 as a real value, so it becomes `nan` here and
      the renderer masks on it. */

      double vcs = interpolate (cs, xp, yp);
      double vf  = interpolate (f, xp, yp);
      double vd  = interpolate (D2c, xp, yp);
      double vu  = interpolate (vel, xp, yp);
      double vax = interpolate (ux, xp, yp);
      fprintf (ferr, "%g %g %g %g %g %g %g\n", xp, yp,
               vcs > 0.5*nodata ? nan("") : vcs,
               vf  > 0.5*nodata ? nan("") : vf,
               vd  > 0.5*nodata ? nan("") : vd,
               vu  > 0.5*nodata ? nan("") : vu,
               fabs(vax) > 0.5*nodata ? nan("") : vax);
    }
  }
  return 0;
}
