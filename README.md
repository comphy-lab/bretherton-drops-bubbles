# Bretherton drops and bubbles

Axisymmetric Basilisk simulations of long drops and bubbles translating in
liquid-filled capillary tubes, aimed at recovering the classical Bretherton
film-thickness law and then comparing drop and bubble counterparts.

## Overview

A long gas bubble or immiscible drop pushed through a liquid-filled tube
never touches the wall: it rides on a thin film of the continuous phase
deposited by its front meniscus. For a bubble at small capillary number
$\mathrm{Ca}_b = \mu_c U_b/\sigma$, lubrication theory
([Bretherton, 1961](https://doi.org/10.1017/S0022112061000160)) gives

$$\frac{b}{R_{tube}} \simeq 1.34\,\mathrm{Ca}_b^{2/3}, \qquad
W = \frac{U_b - U}{U_b} \simeq 1.29\,(3\,\mathrm{Ca}_b)^{2/3},$$

where $b$ is the film thickness, $U$ the mean speed of the carrying
liquid and $U_b$ the bubble speed. At moderate $\mathrm{Ca}_b$ the
Aussillous & Quéré (2000) fit
$b/R_{tube} = 1.34\,\mathrm{Ca}_b^{2/3}/(1 + 3.35\,\mathrm{Ca}_b^{2/3})$
describes Taylor's data.

The tube wall is an **embedded boundary** (`embed.h`), fully wetted by
the continuous phase, and the interface is tracked with **VOF**
(`two-phase.h`). The solver is the axisymmetric incompressible
Navier–Stokes solver (`axi.h` + `navier-stokes/centered.h`) with surface
tension (`tension.h`) and momentum-conserving VOF advection
(`navier-stokes/conserving.h`).

## Non-dimensionalisation

Repeating variables: surface tension $\sigma$, continuous-phase dynamic
viscosity $\mu_c$ and the volume-equivalent drop/bubble radius
$R = (3V/4\pi)^{1/3}$. Hence lengths are in units of $R$, velocities in
units of the visco-capillary velocity $\sigma/\mu_c$, time in units of
$\mu_c R/\sigma$ and pressure in units of $\sigma/R$. A dimensionless
velocity **is** a capillary number, so the measured tip velocity of the
drop/bubble is directly $\mathrm{Ca}_b$.

Runtime control parameters (see `default.params`):

| Key | Meaning | Bubble default |
|-----|---------|----------------|
| `Ca` | imposed mean inlet capillary number $\mu_c U/\sigma$ | 0.05 |
| `La` | Laplace number $\rho_c\sigma R/\mu_c^2$ (keep $La\,Ca \ll 1$ for the visco-capillary regime) | 1 |
| `muR` | viscosity ratio $\mu_d/\mu_c$ | 0.01 |
| `rhoR` | density ratio $\rho_d/\rho_c$ | 0.001 |
| `Rtube` | tube radius in units of $R$ ($<1$ confines the drop/bubble) | 0.7 |
| `MAXlevel` | finest grid level (the film needs several cells: $b \gtrsim 4\,L_{domain}/2^{MAXlevel}$) | 10 |

## Embedded boundaries + VOF: the incompatibility and its handling

`embed.h` and `two-phase.h` are not fully compatible out of the box.
The couplings and their treatment (all collected in
`src-local/embed-vof-tube.h`) are:

1. **Stale axi+embed metric** (upstream test `src/test/missing_metric.c`):
   with `AXI` and `EMBED` the metric is $c_m = y\,c_s$, $f_m = y\,f_s$,
   but `axi.h` computes it once at startup. Every `solid()` call and
   every `adapt_wavelet()` must be followed by
   `cm_update()`/`fm_update()`/`restriction()` —
   wrapped here as `embed_axi_metric_sync()`.
2. **VOF fraction leaking into the solid**: grid adaptation prolongates
   `f` without knowledge of `cs`, and the height-function curvature
   used by `tension.h` is not embed-aware, so spurious `f` inside the
   wall corrupts the film curvature. `vof_solid_cleanup()` resets `f`
   to the continuous-phase value in full-solid cells after adaptation.
3. **Interface–wall separation**: `heights.h`/`curvature.h` ignore
   `cs`, so the wetting film must stay resolved by several cells; the
   case adapts on `cs` to keep the wall region refined and logs a
   warning if the film thins below one fine cell.

Surface tension itself is embed-safe (`iforce.h` skips faces with
$f_m = 0$), and `vof.h` carries explicit `EMBED` branches for the
advection. Both couplings are verified in `testCases/`.

## Requirements

- [Basilisk C](http://basilisk.fr) (`qcc` in `PATH`), e.g. from
  [comphy-lab/basilisk-C](https://github.com/comphy-lab/basilisk-C)
- `bash`, `awk`, `python3` (post-processing)
- OpenMP optional (`--threads N`)

## Quick start

```bash
# smoke + verification tests (~1 minute on a laptop)
bash testCases/runSmokeTests.sh

# single case with the defaults (bubble, Ca = 0.05)
bash runSimulation.sh

# single case, custom parameter file, 4 OpenMP threads
bash runSimulation.sh myCase.params --threads 4

# bubble validation sweep (5 cases across Ca; production hardware)
bash runParameterSweep.sh --dry-run     # inspect first
bash runParameterSweep.sh --threads 4

# drop counterpart sweep
bash runParameterSweep.sh sweep-drop.params --threads 4

# film thickness and Ca_b against Bretherton / Aussillous-Quere
python3 postProcess/bretherton_film.py simulationCases/10?? --out film.csv
```

## Repository structure

```
├── simulationCases/bretherton.c - axisymmetric drop/bubble in an embedded tube (main case)
├── src-local/embed-vof-tube.h - embed + axi + VOF compatibility layer
├── src-local/params.h - typed runtime-parameter accessors with defaults
├── src-local/parse_params.h - low-level key=value parameter parser
├── testCases/embedAxiVofAdvection.c - verification: embed+axi+VOF advection vs exact solution
├── testCases/laplaceEmbedTube.c - verification: Young-Laplace balance inside an embedded tube
├── testCases/runSmokeTests.sh - compile-and-run smoke/verification driver
├── testCases/smoke.params - coarse short-run parameters for the smoke test
├── testCases/README.md - test classification (smoke / verification / validation)
├── postProcess/getFacets.c - extract interface facets from a snapshot
├── postProcess/bretherton_film.py - film thickness and Ca_b from case logs
├── default.params - baseline runtime parameters (bubble)
├── sweep.params - bubble validation sweep over Ca
├── sweep-drop.params - drop counterpart sweep over Ca
├── runSimulation.sh - compile and run one case in simulationCases/<CaseNo>/
└── runParameterSweep.sh - Cartesian SWEEP_* runner built on runSimulation.sh
```

## Outputs

Each case writes to `simulationCases/<CaseNo>/`:

- `c<CaseNo>-log`: per-step diagnostics
  (`i dt t ke dVol/Vol0 xTipF xTipR bFilm`),
- `intermediate/snapshot-*`: time-stamped dumps,
- `restart`: checkpoint for resuming.

`postProcess/bretherton_film.py` fits the front-tip velocity (which is
$\mathrm{Ca}_b$ in code units) and the steady film thickness over the
trailing time window and compares them with the Bretherton and
Aussillous–Quéré predictions.
