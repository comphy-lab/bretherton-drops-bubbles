#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "../src-local/central-film.h"

static int failures;

#define CHECK(condition) do {                                                \
  if (!(condition)) {                                                       \
    fprintf (stderr, "FAIL %s:%d: %s\n", __func__, __LINE__, #condition);  \
    failures++;                                                             \
  }                                                                         \
} while (0)

static CentralFilmConfig standard_config (void)
{
  return (CentralFilmConfig) {
    .advance_window = 1.,
    .minimum_observed_advance = 3.,
    .film_relative_tolerance = 0.02,
    .front_speed_relative_tolerance = 0.02,
    .minimum_coverage = 0.8,
    .maximum_spatial_spread = 0.1,
    .minimum_fresh_fraction = 0.75,
    .minimum_cells = 4.,
    .shape_relative_tolerance = 0.02,
    .hold_windows = 2
  };
}

static CentralFilmMeasurement measurement (double time, double front,
                                           double film)
{
  return (CentralFilmMeasurement) {
    .time = time,
    .front = front,
    .rear = front - 2.,
    .centroid = front - 1.,
    .film = film,
    .spatial_spread = 0.01,
    .coverage = 0.95,
    .fresh_fraction = 0.9,
    .minimum_cells = 6.
  };
}

static bool push (CentralFilmObserver * observer, CentralFilmWindow * window,
                  double time, double front, double film)
{
  CentralFilmMeasurement sample = measurement (time, front, film);
  return central_film_observer_push (observer, sample, window);
}

static void test_quantiles (void)
{
  double odd[] = {9., 1., 5.};
  double even[] = {4., 1., 3., 2.};
  double quartiles[] = {0., 10., 20., 30., 40.};
  CHECK (central_film_median (odd, 3) == 5.);
  CHECK (central_film_median (even, 4) == 2.5);
  CHECK (central_film_quantile (quartiles, 5, 0.25) == 10.);
  CHECK (isnan (central_film_quantile (NULL, 3, 0.5)));
  CHECK (isnan (central_film_quantile (quartiles, 0, 0.5)));
  CHECK (isnan (central_film_quantile (quartiles, 5, -0.01)));
  CHECK (isnan (central_film_quantile (quartiles, 5, 1.01)));
  CHECK (isnan (central_film_quantile (quartiles, 5, NAN)));
  CHECK (isnan (central_film_quantile (quartiles, 5, INFINITY)));
}

static void test_stable_translating_film_converges (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());

  CHECK (!push (&observer, &window, 0., 0., 0.1));
  CHECK (push (&observer, &window, 1., 1., 0.1));
  CHECK (window.quality_ok);
  CHECK (!window.stable);
  CHECK (!window.film_converged);
  CHECK (push (&observer, &window, 2., 2., 0.1));
  CHECK (!window.stable);
  CHECK (window.hold_count == 0);
  CHECK (!window.film_converged);
  CHECK (push (&observer, &window, 3., 3., 0.1));
  CHECK (window.stable);
  CHECK (window.hold_count == 1);
  CHECK (!window.film_converged);
  CHECK (push (&observer, &window, 4., 4., 0.1));
  CHECK (window.stable);
  CHECK (window.hold_count == 2);
  CHECK (window.film_converged);
  CHECK (window.shape_steady);
}

static void test_reinitialisation_uses_relative_advance (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());

  CHECK (!push (&observer, &window, 1000., 10000., 0.1));
  CHECK (!window.film_converged);
  CHECK (push (&observer, &window, 1001., 10001., 0.1));
  CHECK (!window.stable);
  CHECK (!window.film_converged);
  CHECK (push (&observer, &window, 1002., 10002., 0.1));
  CHECK (!window.film_converged);
  CHECK (push (&observer, &window, 1003., 10003., 0.1));
  CHECK (!window.film_converged);
  CHECK (push (&observer, &window, 1004., 10004., 0.1));
  CHECK (window.film_converged);
}

static void check_bad_quality_sample (CentralFilmMeasurement bad)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());
  CHECK (!central_film_observer_push (&observer, bad, &window));
  CHECK (!observer.initialized);
  CHECK (!observer.have_previous_window);
  CHECK (observer.hold_count == 0);
}

static void test_quality_screens (void)
{
  CentralFilmMeasurement bad = measurement (0., 0., 0.1);
  bad.fresh_fraction = 0.2;
  check_bad_quality_sample (bad);

  bad = measurement (0., 0., 0.1);
  bad.coverage = 0.5;
  check_bad_quality_sample (bad);

  bad = measurement (0., 0., 0.1);
  bad.minimum_cells = 2.;
  check_bad_quality_sample (bad);

  bad = measurement (0., 0., 0.1);
  bad.spatial_spread = 0.2;
  check_bad_quality_sample (bad);

  bad = measurement (0., 0., 0.);
  check_bad_quality_sample (bad);
}

static void test_bad_quality_breaks_a_hold (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  CentralFilmConfig config = standard_config();
  CentralFilmMeasurement sample;
  central_film_observer_init (&observer, config);

  CHECK (!push (&observer, &window, 0., 0., 0.1));
  CHECK (push (&observer, &window, 1., 1., 0.1));
  CHECK (push (&observer, &window, 2., 2., 0.1));
  CHECK (window.hold_count == 0);
  CHECK (push (&observer, &window, 3., 3., 0.1));
  CHECK (window.hold_count == 1);

  sample = measurement (4., 4., 0.1);
  sample.coverage = 0.5;
  CHECK (!central_film_observer_push (&observer, sample, &window));
  CHECK (!observer.initialized);
  CHECK (!observer.have_previous_window);
  CHECK (observer.hold_count == 0);
  CHECK (!push (&observer, &window, 5., 5., 0.1));
  CHECK (push (&observer, &window, 6., 6., 0.1));
  CHECK (!window.stable);
  CHECK (!window.film_converged);
}

static void test_drifting_film_is_not_stable (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());
  CHECK (!push (&observer, &window, 0., 0., 0.1));
  CHECK (push (&observer, &window, 1., 1., 0.1));
  CHECK (push (&observer, &window, 2., 2., 0.2));
  CHECK (window.film_drift > observer.config.film_relative_tolerance);
  CHECK (!window.stable);
  CHECK (!window.film_converged);
}

static void test_changing_front_speed_is_not_stable (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());
  CHECK (!push (&observer, &window, 0., 0., 0.1));
  CHECK (push (&observer, &window, 1., 1., 0.1));
  CHECK (push (&observer, &window, 3., 2., 0.1));
  CHECK (window.front_speed_drift >
         observer.config.front_speed_relative_tolerance);
  CHECK (!window.stable);
  CHECK (!window.film_converged);
}

static void test_film_and_shape_milestones_are_distinct (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());

  for (int i = 0; i <= 4; i++) {
    CentralFilmMeasurement sample = measurement ((double) i, (double) i, 0.1);
    sample.rear = -2. + 0.5*i;
    sample.centroid = -1. + 0.8*i;
    bool completed = central_film_observer_push (&observer, sample, &window);
    CHECK (completed == (i > 0));
  }
  CHECK (window.film_converged);
  CHECK (!window.shape_steady);
  CHECK (!window.shape_converged);
}

static void test_shape_convergence_needs_consecutive_windows (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());

  for (int i = 0; i <= 4; i++) {
    CentralFilmMeasurement sample = measurement ((double) i, (double) i, 0.1);
    sample.rear = -2. + 0.5*i;
    sample.centroid = -1. + 0.8*i;
    (void) central_film_observer_push (&observer, sample, &window);
  }
  CHECK (window.film_converged);
  CHECK (!window.shape_steady);

  CentralFilmMeasurement steady = measurement (5., 5., 0.1);
  steady.rear = 1.;
  steady.centroid = 3.2;
  CHECK (central_film_observer_push (&observer, steady, &window));
  CHECK (window.film_converged);
  CHECK (window.shape_steady);
  CHECK (window.shape_hold_count == 1);
  CHECK (!window.shape_converged);

  steady = measurement (6., 6., 0.1);
  steady.rear = 2.;
  steady.centroid = 4.2;
  CHECK (central_film_observer_push (&observer, steady, &window));
  CHECK (window.shape_steady);
  CHECK (window.shape_hold_count == 2);
  CHECK (window.shape_converged);
}

static void assert_discontinuity_restarts_observation
  (CentralFilmMeasurement discontinuity, bool discontinuity_is_new_anchor)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  central_film_observer_init (&observer, standard_config());
  CHECK (!push (&observer, &window, 0., 0., 0.1));
  CHECK (push (&observer, &window, 1., 1., 0.1));
  CHECK (push (&observer, &window, 2., 2., 0.1));
  CHECK (window.hold_count == 0);
  CHECK (push (&observer, &window, 3., 3., 0.1));
  CHECK (window.hold_count == 1);

  CHECK (!central_film_observer_push (&observer, discontinuity, &window));

  /* A bad observation breaks continuity.  The next valid sample is a new
     anchor and cannot close a window using history from before the gap. */
  CHECK (push (&observer, &window, 20., 20., 0.1) ==
         discontinuity_is_new_anchor);
  if (!discontinuity_is_new_anchor)
    CHECK (push (&observer, &window, 21., 21., 0.1));
  CHECK (!window.stable);
  CHECK (window.hold_count == 0);
  CHECK (!window.film_converged);
}

static void test_nonfinite_sample_breaks_history (void)
{
  CentralFilmMeasurement bad = measurement (3., 3., 0.1);
  bad.film = NAN;
  assert_discontinuity_restarts_observation (bad, false);

  bad = measurement (3., 3., 0.1);
  bad.front = INFINITY;
  assert_discontinuity_restarts_observation (bad, false);
}

static void test_nonmonotonic_time_breaks_history (void)
{
  CentralFilmMeasurement bad = measurement (1.5, 3., 0.1);
  assert_discontinuity_restarts_observation (bad, true);

  bad = measurement (3., 4., 0.1);
  assert_discontinuity_restarts_observation (bad, true);
}

static void test_nonmonotonic_front_breaks_history (void)
{
  CentralFilmMeasurement bad = measurement (3., 1.5, 0.1);
  assert_discontinuity_restarts_observation (bad, true);

  bad = measurement (4., 3., 0.1);
  assert_discontinuity_restarts_observation (bad, true);
}

static void test_null_arguments_are_rejected (void)
{
  CentralFilmObserver observer;
  CentralFilmWindow window = {0};
  CentralFilmMeasurement sample = measurement (0., 0., 0.1);
  central_film_observer_init (&observer, standard_config());
  CHECK (!central_film_observer_push (NULL, sample, &window));
  CHECK (!central_film_observer_push (&observer, sample, NULL));
}

int main (void)
{
  test_quantiles();
  test_stable_translating_film_converges();
  test_reinitialisation_uses_relative_advance();
  test_quality_screens();
  test_bad_quality_breaks_a_hold();
  test_drifting_film_is_not_stable();
  test_changing_front_speed_is_not_stable();
  test_film_and_shape_milestones_are_distinct();
  test_shape_convergence_needs_consecutive_windows();
  test_nonfinite_sample_breaks_history();
  test_nonmonotonic_time_breaks_history();
  test_nonmonotonic_front_breaks_history();
  test_null_arguments_are_rejected();

  if (failures) {
    fprintf (stderr, "central-film observer: %d assertion(s) failed\n", failures);
    return 1;
  }
  puts ("central-film observer: PASS");
  return 0;
}
