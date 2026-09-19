# Execution checks

This directory holds observer unit tests, a serial smoke test and a serial/MPI
integration check.

The serial smoke test checks that `simulationCases/bretherton.c` compiles and
integrates a few dozen time steps at coarse resolution without producing
a truncated log or blowing up. It compares against **nothing**, so it
supports no claim about the discretisation or about the physics. Its
value is that it fails loudly and in seconds when a change breaks the
build or the parameter path.

| File | Role |
|------|------|
| `smoke.params` | coarse, short-run parameters (`MAXlevel = 9`, `tmax = 0.05`) |
| `runSmokeTests.sh` | builds and runs the case, then checks the log |
| `runParallelSmokeTests.py` | compares serial and two-rank MPI runs and restarts |
| `centralFilmObserver.c` | tests convergence windows, quality gates and history resets |
| `test_central_film_observer.py` | same observer contract in Python, plus log stitching across a restart |
| `runCentralFilmTests.sh` | compiles the C tests and runs the Python log-acceptance tests |
| `regridTube.c` | measures transfer and projection changes when refining a saved tube state |

Cases that compare against an exact solution live in
[`../verificationCases/`](../verificationCases/) and are driven by
`../verificationCases/runVerification.sh`. The campaign that compares
against independent experimental data — Taylor (1961) and
Aussillous & Quéré (2000) — is configured by `sweep.params` and
`sweep-drop.params` at the repository root and is meant for production
hardware. Those three evidence sources support three different claims
and are deliberately kept apart.

## Running

```bash
bash testCases/runSmokeTests.sh   # this directory only
bash runTests.sh --unit          # C observer tests and Python log stitching, without Basilisk
bash runTests.sh                 # observer, verification and smoke tests
```

The build goes to `simulationCases/9999/`, which is gitignored.

The parallel integration check needs Linux with Python pidfd support,
`qcc`, `mpicc`, and an Open MPI rankfile naming two allocated cores.
It builds serial and MPI modes once, compares short fresh and restarted runs,
checks their logs and snapshots, and exercises collective early termination.
It also compares refined restarts and checks that rejecting a changed domain
length leaves the seed checkpoint intact. Short runs must report the expected
incomplete terminal condition; a nonzero exit alone does not pass.
Run it only inside the corresponding reserved CPU set with a healthy MPI
toolchain:

```bash
python3 testCases/runParallelSmokeTests.py --rankfile /absolute/path/to/rankfile \
  --work-root /absolute/path/to/empty-test-directory
```

The check retains command logs and stops on failed execution, VOF CFL warnings
or disagreement beyond its stated tolerances. It checks parallel execution
and restart consistency; performance and mesh convergence require separate
tests.
