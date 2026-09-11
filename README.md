# Bretherton drops and bubbles

Axisymmetric Basilisk solver for a gas bubble or immiscible drop translating
through a liquid-filled circular tube. Runtime parameter files control the
case; the current validation programme concerns Newtonian bubbles.

The scientific report is in `docs/Newtonian-Validation/`: `main.tex`,
`references.bib`, figures, compact comparison data and the compiled `main.pdf`.

## Requirements

- A project-local [Basilisk C](https://github.com/comphy-lab/basilisk-C)
  installation; the validation campaign uses `v2026-08-30`.
- Bash, a C compiler and `qcc`. The runners load `.project_config` when present.
- OpenMP for `--threads N`; Python with NumPy for flat-film post-processing.
- Open MPI (`mpicc` and `mpirun`) for `--ranks N`.

## Run a case or sweep

```bash
# Compile and run the coarse smoke test.
bash runTests.sh --smoke

# Run a single case with the bubble defaults.
bash runSimulation.sh default.params --threads 4

# MPI: the rankfile must map ranks to the allocated physical cores.
bash runSimulation.sh default.params --ranks 4 --rankfile /path/to/rankfile

# Inspect the eight-case grid study, then run it.
bash runParameterSweep.sh grid-sensitivity.params --dry-run
bash runParameterSweep.sh grid-sensitivity.params --threads 4 --parallel 2
```

[`default.params`](default.params) defines the physical inputs and numerical
controls. `Ca` is the imposed mean inlet capillary number; the bubble speed is
measured from the output. [`grid-sensitivity.params`](grid-sensitivity.params)
repeats the four inlet values at maximum levels 9 and 11. The separate
[`grid-sensitivity-fine.params`](grid-sensitivity-fine.params) specifies the
level-12 follow-up.

Case output defaults to `simulationCases/<CaseNo>/`. Set `OUTPUT_ROOT` to a
separate run directory for production calculations. Each case contains its
parameters, source, executable, diagnostic log, `restart` dump and
`intermediate/snapshot-*` files. An existing restart resumes that case.
Success requires stable deposited central-film measurements and, by default,
agreement of the front, rear and centroid speeds. Reaching `tmax` or the
outlet buffer returns an explicit incomplete status with a nonzero exit code.
For a refined restart, use a separate case directory containing a copy of the
seed dump, set `freshFront` to its front position, and specify `regridBurnR`.
Changing `Ldomain` while restoring a dump is rejected; start a fresh case to
change the domain length.

For a single-node MPI allocation, `--pe-list 0,1,2,3` can replace the rankfile
when those are the allocated Open MPI core indices. Set `MPIEXEC` to select
another launcher executable. `--build-only` compiles without running;
`--no-build` reuses an executable after checking its source, local headers,
Basilisk lock and binary checksums. `--openmp --threads 1` retains an OpenMP
build for one-thread tests.

## Analysis and numerical checks

[`postProcess/bretherton_flat_film.py`](postProcess/bretherton_flat_film.py)
extracts the flat film and bubble speed from saved interfaces; it uses the
compiled [`getFacets.c`](postProcess/getFacets.c) helper.
[`postProcess/bretherton_film.py`](postProcess/bretherton_film.py) instead reads
the minimum film reported in the case log. The validation report explains the
difference between these measurements.

```bash
# Exact-solution checks followed by the smoke test.
bash runTests.sh
```

The [verification cases](verificationCases/README.md) state their comparators
and tolerances. The [smoke test](testCases/README.md) checks compilation and
short-time execution.

## Source layout

```
├── simulationCases/bretherton.c - production solver
├── src-local/embed-vof-tube.h - embedded-wall and VOF compatibility functions
├── src-local/params.h - typed runtime parameters
├── postProcess/ - interface extraction, film analysis and videos
├── verificationCases/ - exact-solution numerical checks
├── testCases/ - smoke test
└── docs/Newtonian-Validation/ - LaTeX report, bibliography, figures and PDF
```

## Compile the scientific report

Requires `latexmk`, a LaTeX distribution and BibTeX:

```bash
make -C docs/Newtonian-Validation
```

This produces `docs/Newtonian-Validation/main.pdf`. The folder's `README.md`
lists its sources and the optional figure-regeneration command.

## Build and preview the code documentation

```bash
bash .github/scripts/build.sh
bash .github/scripts/deploy.sh
```

The build generates the site in `.github/docs/`; `deploy.sh` starts a local
preview. See [the documentation workflow](.github/Website-generator-readme.md)
for source discovery and GitHub Pages publication.
