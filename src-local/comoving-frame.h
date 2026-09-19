/**
# Co-moving frame controller and renewal-window observer

Solver-independent helpers for a window that translates with the bubble.

## Frame transformation

With lab velocity $\mathbf{v} = \mathbf{u} + U(t)\,\mathbf{e}_x$, the
equations for $\mathbf{u}$ are the Navier--Stokes equations plus one uniform
body force per unit mass, $-\dot U\,\mathbf{e}_x$, acting on both phases.
The transformation is exact, so $U(t)$ cannot alter the laboratory physics;
it only decides where the bubble sits in the window. The caller applies
`-dUdt` through the Basilisk acceleration face vector and shifts every
velocity boundary condition by `-U`.

## Closure

Basilisk has no global unknown, so $U$ is set by a controller on the
dispersed-phase centroid $x_c$. With $e = x_c - x_{target}$ and
$\dot x_c = U_b - U$,
$$
\dot U = \frac{2}{\tau}\,\dot x_c + \frac{1}{\tau^2}\,e
\qquad\Rightarrow\qquad
\ddot e + \frac{2}{\tau}\dot e + \frac{e}{\tau^2} = \dot U_b ,
$$
a critically damped loop forced only by the rate of change of the true
bubble speed. Once $U_b$ is steady, $e \to 0$ and $U \to U_b$. Galilean
invariance decouples the loop from the physics; $\tau$ shapes the transient
and the size of the fictitious force only.

The controller settles the impulse of the previous step before it computes
the next acceleration: the momentum received $\dot U\,\Delta t$ during the
step just taken, so the frame speed used by the boundary conditions advances
by the same amount. At steady state $\dot U = 0$ and the two are exactly
consistent.

## Renewal-window observer

Steady state is judged in film renewals $L_b/U$ because the film is
material: the central 30--70 % window is fresh only after the old film has
been swept to the rear. Successive windows must agree on the film, the frame
speed, the centroid offset and the bubble length for `hold_windows` windows.
*/

#ifndef COMOVING_FRAME_H
#define COMOVING_FRAME_H

#include <math.h>
#include <stdbool.h>
#include <string.h>

typedef struct {
  double tau;            // controller time constant (visco-capillary units)
  double target;         // centroid target in the window
  bool prescribed;       // true: U follows prescribed_U, controller off
  double prescribed_U;
  double U, dUdt, v_rel, error;
  double displacement;   // integral of U since t_change (frame travel)
  double xc_prev, t_prev;
  bool have_prev;
} ComovingFrame;

static inline void comoving_frame_init (ComovingFrame * fr, double U0,
                                        double target, double tau,
                                        bool prescribed)
{
  memset (fr, 0, sizeof(*fr));
  fr->tau = tau;
  fr->target = target;
  fr->prescribed = prescribed;
  fr->prescribed_U = U0;
  fr->U = U0;
}

/**
Advances the frame with the measured centroid at time `t`. Returns the frame
acceleration to impose over the coming step. The first call only records the
state. A non-advancing clock resets the history rather than dividing by zero.
*/
static inline double comoving_frame_update (ComovingFrame * fr, double xc,
                                            double t)
{
  if (!fr->have_prev) {
    fr->have_prev = true;
    fr->xc_prev = xc;
    fr->t_prev = t;
    fr->error = xc - fr->target;
    fr->v_rel = 0.;
    fr->dUdt = 0.;
    return 0.;
  }
  double dt = t - fr->t_prev;
  if (!(dt > 0.)) {
    fr->xc_prev = xc;
    fr->t_prev = t;
    return fr->dUdt;
  }
  // settle the impulse of the step just taken
  fr->U += fr->dUdt*dt;
  fr->displacement += fr->U*dt;
  fr->v_rel = (xc - fr->xc_prev)/dt;
  fr->error = xc - fr->target;
  if (fr->prescribed) {
    fr->dUdt = 0.;
    fr->U = fr->prescribed_U;
  }
  else
    fr->dUdt = (2./fr->tau)*fr->v_rel + fr->error/(fr->tau*fr->tau);
  fr->xc_prev = xc;
  fr->t_prev = t;
  return fr->dUdt;
}

/** Windowed stationarity in renewal units. */

typedef struct {
  double window_renewals, burn_renewals;
  double film_tol, speed_tol, position_tol, length_tol;
  double minimum_coverage, maximum_spread, minimum_fresh, minimum_cells;
  int hold_windows;
} RenewalConfig;

typedef struct {
  double time, film, U, xc, length, error;
  double spread, coverage, fresh, cells;
} RenewalSample;

typedef struct {
  double t_start, t_end, duration;
  double film, U, xc, length, error_abs_max;
  double film_drift, speed_drift, length_drift;
  double min_coverage, max_spread, min_fresh, min_cells;
  int hold_count;
  bool quality_ok, stable, converged;
} RenewalWindow;

typedef struct {
  RenewalConfig config;
  bool initialized, have_previous, burn_done;
  double t_burn_start, t_window_start, prev_time;
  double film_int, U_int, xc_int, len_int, elapsed;
  double err_max, min_cov, max_spr, min_fresh, min_cells;
  RenewalSample prev;
  double prev_film, prev_U, prev_len;
  int hold_count;
} RenewalObserver;

static inline void renewal_observer_init (RenewalObserver * ob,
                                          RenewalConfig c)
{
  memset (ob, 0, sizeof(*ob));
  ob->config = c;
}

static inline void renewal_window_reset (RenewalObserver * ob,
                                         RenewalSample s)
{
  ob->t_window_start = s.time;
  ob->prev = s;
  ob->film_int = ob->U_int = ob->xc_int = ob->len_int = ob->elapsed = 0.;
  ob->err_max = fabs(s.error);
  ob->min_cov = s.coverage; ob->max_spr = s.spread;
  ob->min_fresh = s.fresh; ob->min_cells = s.cells;
}

/**
Pushes one sample. `renewal` is the current estimate of $L_b/U$. Returns true
when a window has completed and fills `w`. Samples failing the quality gates
reset the current window but not the burn clock.
*/
static inline bool renewal_observer_push (RenewalObserver * ob,
                                          RenewalSample s, double renewal,
                                          RenewalWindow * w)
{
  const RenewalConfig * c = &ob->config;
  if (!isfinite(s.time) || !isfinite(s.film) || !isfinite(s.U) ||
      !isfinite(s.xc) || !isfinite(s.length) || !(renewal > 0.))
    return false;
  if (!ob->initialized) {
    ob->initialized = true;
    ob->t_burn_start = s.time;
    renewal_window_reset (ob, s);
    return false;
  }
  if (!ob->burn_done) {
    if (s.time - ob->t_burn_start < c->burn_renewals*renewal) {
      ob->prev = s;
      return false;
    }
    ob->burn_done = true;
    renewal_window_reset (ob, s);
    return false;
  }
  bool quality = s.film > 0. && s.coverage >= c->minimum_coverage &&
    s.spread <= c->maximum_spread && s.fresh >= c->minimum_fresh &&
    (c->minimum_cells <= 0. || s.cells >= c->minimum_cells);
  if (!quality) {
    renewal_window_reset (ob, s);
    ob->hold_count = 0;
    return false;
  }
  double dt = s.time - ob->prev.time;
  if (dt > 0.) {
    ob->film_int += .5*(s.film + ob->prev.film)*dt;
    ob->U_int += .5*(s.U + ob->prev.U)*dt;
    ob->xc_int += .5*(s.xc + ob->prev.xc)*dt;
    ob->len_int += .5*(s.length + ob->prev.length)*dt;
    ob->elapsed += dt;
  }
  ob->err_max = fmax(ob->err_max, fabs(s.error));
  ob->min_cov = fmin(ob->min_cov, s.coverage);
  ob->max_spr = fmax(ob->max_spr, s.spread);
  ob->min_fresh = fmin(ob->min_fresh, s.fresh);
  ob->min_cells = fmin(ob->min_cells, s.cells);
  ob->prev = s;
  if (s.time - ob->t_window_start < c->window_renewals*renewal ||
      ob->elapsed <= 0.)
    return false;

  memset (w, 0, sizeof(*w));
  w->t_start = ob->t_window_start; w->t_end = s.time;
  w->duration = s.time - ob->t_window_start;
  w->film = ob->film_int/ob->elapsed;
  w->U = ob->U_int/ob->elapsed;
  w->xc = ob->xc_int/ob->elapsed;
  w->length = ob->len_int/ob->elapsed;
  w->error_abs_max = ob->err_max;
  w->min_coverage = ob->min_cov; w->max_spread = ob->max_spr;
  w->min_fresh = ob->min_fresh; w->min_cells = ob->min_cells;
  w->quality_ok = w->film > 0. && w->U > 0.;
  if (ob->have_previous) {
    w->film_drift = fabs(w->film - ob->prev_film)/fabs(ob->prev_film);
    w->speed_drift = fabs(w->U - ob->prev_U)/fabs(ob->prev_U);
    w->length_drift = fabs(w->length - ob->prev_len)/fabs(ob->prev_len);
    w->stable = w->quality_ok &&
      w->film_drift <= c->film_tol && w->speed_drift <= c->speed_tol &&
      w->length_drift <= c->length_tol &&
      w->error_abs_max <= c->position_tol;
    ob->hold_count = w->stable ? ob->hold_count + 1 : 0;
  }
  else {
    w->film_drift = w->speed_drift = w->length_drift = INFINITY;
    ob->have_previous = true;
    ob->hold_count = 0;
  }
  ob->prev_film = w->film; ob->prev_U = w->U; ob->prev_len = w->length;
  w->hold_count = ob->hold_count;
  w->converged = ob->hold_count >= c->hold_windows;
  renewal_window_reset (ob, s);
  return true;
}

#endif
