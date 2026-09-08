/**
# embed-vof-tube.h

Compatibility layer for combining embedded boundaries (`embed.h`) with
axisymmetric Volume-Of-Fluid two-phase flow (`axi.h` + `two-phase.h`).

## Why this header exists

Basilisk's embedded boundaries and VOF advection are *almost* compatible
out of the box, but three couplings must be handled by the case:

1. **Missing axi metric after `solid()` and after adaptation.** With both
   `AXI` and `EMBED` defined, the metric fields must satisfy
   $c_m = y\,c_s$ and $f_m = y\,f_s$. `axi.h` computes them once in its
   `metric` event (when `cs = 1` everywhere), so every later call to
   `solid()` — and every `adapt_wavelet()` on trees, which re-prolongates
   `cs`/`fs` — leaves `cm`/`fm` stale. `axi.h` provides `cm_update()` and
   `fm_update()` for exactly this purpose and the upstream test
   [`src/test/missing_metric.c`](http://basilisk.fr/src/test/missing_metric.c)
   guards the coupling. `embed_axi_metric_sync()` wraps the required calls.

2. **VOF fraction inside the solid.** `vof.h` only updates cells with
   `cs > 0`, but grid adaptation prolongates `f` without knowledge of
   `cs`, so refined solid cells can inherit non-zero `f` from fluid
   neighbours. The height-function curvature used by `tension.h` is not
   embed-aware: any spurious `f` inside the solid corrupts height columns
   near the wall and hence the film curvature. `vof_solid_cleanup()`
   resets `f` to the continuous-phase value (here `f = 0`) in all
   full-solid cells. This is consistent with a wall perfectly wetted by
   the continuous phase.

3. **Interface–wall separation.** `heights.h`/`curvature.h` have no
   knowledge of `cs`, so the wetting film between the interface and the
   embedded wall must stay resolved by several cells. This header cannot
   enforce that; cases must refine the wall region (e.g. by adapting on
   `cs`) and monitor the film thickness.

Surface tension itself is safe: `iforce.h` guards its face loop with
`fm.x[] > 0`, so no force is applied on embedded faces.

## Public API

- `embed_axi_metric_sync()`: recompute `cm`/`fm` and restrict
  `{cs, fs, cm, fm}`; call after `solid()` and after `adapt_wavelet()`.
- `tube_solid(Rtube)`: embed a cylindrical tube wall of radius `Rtube`
  (fluid at `y < Rtube`) and synchronise the metric.
- `vof_solid_cleanup(f)`: clamp `f` and reset it to the continuous phase
  in full-solid cells; call after `adapt_wavelet()`.
*/

#ifndef EMBED_VOF_TUBE_H
#define EMBED_VOF_TUBE_H

/**
### embed_axi_metric_sync()

Recomputes the metric fields from the current `cs`/`fs` and restricts
all four fields on trees. With `EMBED` but no `AXI`, only the
restriction is needed (the metric is the solid fraction itself, handled
by `embed.h`). After `tube_solid()` registers a fixed cylindrical wall,
the fractions and metrics are reconstructed from that geometry at every
stored tree level, including MPI neighbour cells. This specialisation
assumes the wall remains a straight tube of the registered radius.
*/
#if TREE && AXI && EMBED
static double tube_metric_radius;

// Integral of the radius over the fluid portion of a radial cell interval,
// divided by its width. This is both cm and the axial-face metric of a tube.
static double tube_axial_face_metric (double r, double d)
{
  double lo = r - d/2., hi = min(r + d/2., tube_metric_radius);
  return hi > lo ? (hi*hi - lo*lo)/(2.*d) : 0.;
}
#endif

static inline void embed_axi_metric_sync (void)
{
#if defined(AXI) && defined(EMBED)
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
#if TREE
  // A stationary tube has an exact geometric description at every level.
  // Populate inactive parents and MPI neighbour storage as well as leaves:
  // face restriction can read either side of a partition boundary.
  if (tube_metric_radius > 0.)
    foreach_cell() {
      cs[] = clamp((tube_metric_radius - y)/Delta + .5, 0., 1.);
      fs.x[] = cs[];
      if (allocated(1)) fs.x[1] = cs[];
      fs.y[] = y - Delta/2. < tube_metric_radius ? 1. : 0.;
      if (allocated(0,1)) fs.y[0,1] = y + Delta/2. < tube_metric_radius ? 1. : 0.;
      cm[] = tube_axial_face_metric(y, Delta);
      fm.x[] = tube_axial_face_metric(y, Delta);
      if (allocated(1)) fm.x[1] = tube_axial_face_metric(y, Delta);
      double r = y - Delta/2.;
      fm.y[] = r < tube_metric_radius ? max(r,1e-20) : 0.;
      r = y + Delta/2.;
      if (allocated(0,1)) fm.y[0,1] = r < tube_metric_radius ? max(r,1e-20) : 0.;
    }
#endif
#endif
#if defined(EMBED) && TREE
  // foreach_cell() does not perform automatic boundary-state tracking.
  for (scalar s in {cs, fs, cm, fm}) set_dirty_stencil(s);
  restriction ({cs, fs, cm, fm});
#endif
}

/**
### tube_solid()

Embeds the tube wall: solid for $y > R_{tube}$, fluid below.

#### Parameters
- `Rtube`: tube radius in code units.
*/
static inline void tube_solid (double Rtube)
{
#if TREE && AXI && EMBED
  tube_metric_radius = Rtube;
#endif
  solid (cs, fs, Rtube - y);
  embed_axi_metric_sync();
}

/**
### vof_solid_cleanup()

Enforces the continuous-phase value `f = 0` in full-solid cells and
clamps `f` elsewhere. Cut cells (`0 < cs < 1`) are left untouched: their
volume fraction is a valid fluid quantity.

#### Parameters
- `c`: VOF volume fraction field (dispersed phase at `c = 1`).
*/
static inline void vof_solid_cleanup (scalar c)
{
#if defined(EMBED)
  foreach() {
    if (cs[] <= 0.)
      c[] = 0.;
    else
      c[] = clamp (c[], 0., 1.);
  }
#endif
}

#endif
