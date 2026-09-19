/**
# Central-film observation windows

Small, solver-independent helpers for deciding whether a measured deposited
film is quasi-steady.  The caller owns interface reconstruction and supplies
one spatially screened measurement at a time.  A film gate and a whole-shape
diagnostic are deliberately separate: a steady central film does not imply
that both menisci form a travelling shape.
*/

#ifndef CENTRAL_FILM_H
#define CENTRAL_FILM_H

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double advance_window, minimum_observed_advance;
  double film_relative_tolerance, front_speed_relative_tolerance;
  double minimum_coverage, maximum_spatial_spread;
  double minimum_fresh_fraction, minimum_cells;
  double shape_relative_tolerance;
  int hold_windows;
} CentralFilmConfig;

typedef struct {
  double time, front, rear, centroid;
  double film, spatial_spread, coverage, fresh_fraction, minimum_cells;
} CentralFilmMeasurement;

typedef struct {
  double advance, film, front_speed, rear_speed, centroid_speed;
  double film_drift, front_speed_drift;
  double shape_speed_mismatch, length_drift;
  double minimum_coverage, maximum_spatial_spread;
  double minimum_fresh_fraction, minimum_cells;
  int hold_count, shape_hold_count;
  bool quality_ok, stable, film_converged, shape_steady, shape_converged;
} CentralFilmWindow;

typedef struct {
  CentralFilmConfig config;
  CentralFilmMeasurement first, previous, window_start;
  double film_integral, covered_advance;
  double previous_window_film, previous_front_speed;
  double minimum_coverage, maximum_spatial_spread;
  double minimum_fresh_fraction, minimum_cells;
  int hold_count, shape_hold_count;
  bool initialized, have_previous_window;
} CentralFilmObserver;

static int central_film_compare_double (const void * a, const void * b)
{
  const double da = *(const double *) a, db = *(const double *) b;
  return (da > db) - (da < db);
}

/** Sorts `values` and returns its linearly interpolated quantile. */
static inline double central_film_quantile (double * values, int count,
                                            double quantile)
{
  if (!values || count <= 0 || !isfinite(quantile) ||
      quantile < 0. || quantile > 1.)
    return NAN;
  qsort (values, (size_t) count, sizeof(double),
         central_film_compare_double);
  double index = quantile*(count - 1), lower = floor(index);
  int lo = (int) lower, hi = lo + 1 < count ? lo + 1 : lo;
  return values[lo] + (index - lower)*(values[hi] - values[lo]);
}

static inline double central_film_median (double * values, int count)
{
  return central_film_quantile (values, count, 0.5);
}

static inline void central_film_observer_init (CentralFilmObserver * observer,
                                               CentralFilmConfig config)
{
  memset (observer, 0, sizeof(*observer));
  observer->config = config;
}

/** Clears all observations while preserving the configured predicates. */
static inline void central_film_observer_reset (CentralFilmObserver * observer)
{
  if (!observer)
    return;
  CentralFilmConfig config = observer->config;
  memset (observer, 0, sizeof(*observer));
  observer->config = config;
}

static inline double central_film_relative_change (double current,
                                                   double previous)
{
  return previous != 0. ? fabs(current - previous)/fabs(previous) : INFINITY;
}

static inline void central_film_window_reset (CentralFilmObserver * observer,
                                              CentralFilmMeasurement sample)
{
  observer->window_start = sample;
  observer->previous = sample;
  observer->film_integral = observer->covered_advance = 0.;
  observer->minimum_coverage = sample.coverage;
  observer->maximum_spatial_spread = sample.spatial_spread;
  observer->minimum_fresh_fraction = sample.fresh_fraction;
  observer->minimum_cells = sample.minimum_cells;
}

/**
Adds a measurement.  Returns true only when an advance window has completed;
`window` is unchanged otherwise.  Initialising a new observer anchors all
cumulative distances to the first supplied sample, which makes a restored run
earn fresh observation windows before it can stop.
*/
static inline bool central_film_observer_push (CentralFilmObserver * observer,
                                               CentralFilmMeasurement sample,
                                               CentralFilmWindow * window)
{
  if (!observer || !window)
    return false;
  if (!isfinite(sample.time) ||
      !isfinite(sample.front) || !isfinite(sample.rear) ||
      !isfinite(sample.centroid) || !isfinite(sample.film) ||
      !isfinite(sample.spatial_spread) || !isfinite(sample.coverage) ||
      !isfinite(sample.fresh_fraction) || !isfinite(sample.minimum_cells)) {
    central_film_observer_reset (observer);
    return false;
  }

  const CentralFilmConfig * c = &observer->config;
  bool sampleQuality = sample.film > 0. &&
    sample.coverage >= c->minimum_coverage &&
    sample.spatial_spread <= c->maximum_spatial_spread &&
    sample.fresh_fraction >= c->minimum_fresh_fraction &&
    (c->minimum_cells <= 0. || sample.minimum_cells >= c->minimum_cells);
  if (!sampleQuality) {
    central_film_observer_reset (observer);
    return false;
  }

  if (!observer->initialized) {
    observer->initialized = true;
    observer->first = sample;
    central_film_window_reset (observer, sample);
    return false;
  }

  if (sample.time <= observer->previous.time ||
      sample.front <= observer->previous.front) {
    central_film_observer_reset (observer);
    observer->initialized = true;
    observer->first = sample;
    central_film_window_reset (observer, sample);
    return false;
  }

  double advance = sample.front - observer->previous.front;
  if (advance > 0.) {
    observer->film_integral +=
      0.5*(sample.film + observer->previous.film)*advance;
    observer->covered_advance += advance;
  }
  observer->minimum_coverage = fmin(observer->minimum_coverage,
                                    sample.coverage);
  observer->maximum_spatial_spread = fmax(observer->maximum_spatial_spread,
                                          sample.spatial_spread);
  observer->minimum_fresh_fraction = fmin(observer->minimum_fresh_fraction,
                                          sample.fresh_fraction);
  observer->minimum_cells = fmin(observer->minimum_cells,
                                 sample.minimum_cells);
  observer->previous = sample;

  const double span = sample.front - observer->window_start.front;
  const double elapsed = sample.time - observer->window_start.time;
  if (span < observer->config.advance_window ||
      observer->covered_advance <= 0. || elapsed <= 0.)
    return false;

  memset (window, 0, sizeof(*window));
  window->advance = span;
  window->film = observer->film_integral/observer->covered_advance;
  window->front_speed = span/elapsed;
  window->rear_speed = (sample.rear - observer->window_start.rear)/elapsed;
  window->centroid_speed =
    (sample.centroid - observer->window_start.centroid)/elapsed;
  window->minimum_coverage = observer->minimum_coverage;
  window->maximum_spatial_spread = observer->maximum_spatial_spread;
  window->minimum_fresh_fraction = observer->minimum_fresh_fraction;
  window->minimum_cells = observer->minimum_cells;

  const double length0 = observer->window_start.front -
                         observer->window_start.rear;
  const double length1 = sample.front - sample.rear;
  window->length_drift = length0 > 0. ? fabs(length1 - length0)/length0 :
                                       INFINITY;
  const double speed_scale = fabs(window->centroid_speed);
  window->shape_speed_mismatch = speed_scale > 0. ?
    fmax(fabs(window->front_speed - window->centroid_speed),
         fabs(window->rear_speed - window->centroid_speed))/speed_scale :
    INFINITY;

  window->quality_ok = window->film > 0. && window->front_speed > 0. &&
    window->minimum_coverage >= c->minimum_coverage &&
    window->maximum_spatial_spread <= c->maximum_spatial_spread &&
    window->minimum_fresh_fraction >= c->minimum_fresh_fraction &&
    (c->minimum_cells <= 0. || window->minimum_cells >= c->minimum_cells);
  window->shape_steady =
    window->shape_speed_mismatch <= c->shape_relative_tolerance &&
    window->length_drift <= c->shape_relative_tolerance;

  if (observer->have_previous_window) {
    window->film_drift = central_film_relative_change
      (window->film, observer->previous_window_film);
    window->front_speed_drift = central_film_relative_change
      (window->front_speed, observer->previous_front_speed);
    const double observed = sample.front - observer->first.front;
    window->stable = window->quality_ok &&
      observed >= c->minimum_observed_advance &&
      window->film_drift <= c->film_relative_tolerance &&
      window->front_speed_drift <= c->front_speed_relative_tolerance;
    observer->hold_count = window->stable ? observer->hold_count + 1 : 0;
    observer->shape_hold_count = window->stable && window->shape_steady ?
      observer->shape_hold_count + 1 : 0;
  }
  else {
    window->film_drift = window->front_speed_drift = INFINITY;
    observer->hold_count = 0;
    observer->shape_hold_count = 0;
    observer->have_previous_window = true;
  }

  observer->previous_window_film = window->film;
  observer->previous_front_speed = window->front_speed;
  window->hold_count = observer->hold_count;
  window->shape_hold_count = observer->shape_hold_count;
  window->film_converged = observer->hold_count >= c->hold_windows;
  window->shape_converged = observer->shape_hold_count >= c->hold_windows;
  central_film_window_reset (observer, sample);
  return true;
}

#endif
