# Verification cases

Every case here compares against an **exact solution of the equations the
solver implements**. That supports a claim about the discretisation: that
it converges to the selected mathematical model. It supports no claim
about whether that model describes a real drop or bubble in a tube.

The three evidence sources in this repository are kept apart on purpose:

| Where | Comparator | Claim it supports |
|-------|------------|-------------------|
| `../testCases/` | none — the case only runs | the code compiles and integrates |
| `verificationCases/` (here) | exact solution of the implemented equations | the discretisation converges to the model |
| root `sweep.params`, `sweep-drop.params` | independent experimental data (Taylor 1961; Aussillous & Quéré 2000) | the model reproduces measured behaviour |

Bretherton's $b/R_{tube} = 1.34\,Ca_b^{2/3}$ is an asymptotic solution of
the lubrication limit of the same governing equations, not independent
data. Comparing against it is a limiting-case check whose agreement is
expected to degrade as $Ca_b$ rises; it is not validation.

## Cases

| Case | Comparator | Grid | Cost (serial) |
|------|------------|------|---------------|
| `embedAxiVofAdvection.c` | exact advection $r(t)=\sqrt{r_0^2+2t}$ under $u_r = 1/r$, with an embedded wall | uniform, $N=64$ | seconds |
| `laplaceEmbedTube.c` | Young–Laplace $\Delta p = 2\sigma/R_d$, $\mathbf{u}\equiv 0$; refinement sequence | uniform, $N = 64, 128, 256$ | ~3 min |
| `laplaceEmbedTubeAdapt.c` | same, on an adapted tree, plus direct probes of the `embed-vof-tube.h` couplings | `adapt_wavelet`, levels 4–7 | ~5 s |
| `staticFilmTube.c` | same, with a 4–8 cell interface–wall gap at production property ratios | uniform, $N=256$, 5 runs | ~49 min |

`laplaceEmbedTube.c` and `staticFilmTube.c` write machine-readable raw
errors next to their build: `laplaceEmbedTube-refinement.csv` and
`staticFilmTube-gap.csv`, each with a commented metadata header.

`laplaceEmbedTubeAdapt.c` carries a **negative control**. Built with
`-DSKIP_EMBED_GUARDS` it omits the post-adaptation
`embed_axi_metric_sync()` and `vof_solid_cleanup()` and must fail; the
driver builds that variant automatically and treats a passing guard-less
build as a failure of the suite.

## Running

```bash
bash runTests.sh                          # verification cases, then the smoke test
bash runTests.sh --verification           # this directory only
bash verificationCases/runVerification.sh verificationCases/laplaceEmbedTube.c   # one case

VERIFICATION_THREADS=16 bash runTests.sh  # OpenMP; ~10 min instead of ~49
```

Serial is the default because a recorded verification result should be
reproducible; OpenMP reductions sum in a nondeterministic order and shift
results at the $10^{-5}$ level, which is below every gate here but above
the resolution of some of the comparisons below. Build artefacts go to
`verificationCases/build-*/`, which is gitignored.

## Reference results

Basilisk `qcc`, single core. `embedAxiVofAdvection`, `laplaceEmbedTube`
and `staticFilmTube` run on uniform grids; only `laplaceEmbedTubeAdapt`
adapts. The `staticFilmTube` numbers are from one 16-thread run; see the
reproducibility note below for the spread between independent runs.

`embedAxiVofAdvection` — max relative interface-position error
$2.88\times10^{-5}$ (tolerance $10^{-2}$).

`laplaceEmbedTube` — exact $\Delta p = 5$:

| $N$ | $\Delta$ | $R_d/\Delta$ | rel. $\Delta p$ error | $\max\|\mathbf{u}\|$ |
|---|---|---|---|---|
| 64 | 0.031250 | 12.8 | $2.74\times10^{-3}$ | $1.82\times10^{-4}$ |
| 128 | 0.015625 | 25.6 | $6.05\times10^{-4}$ | $4.14\times10^{-5}$ |
| 256 | 0.0078125 | 51.2 | $1.92\times10^{-4}$ | $1.29\times10^{-5}$ |

Observed order of the pressure-jump error: $p = 2.18$ then $1.66$. Three
levels give two estimates — enough to see a trend, not enough to claim an
asymptotic rate.

`laplaceEmbedTubeAdapt` — P1 (metric consistency) and P2 (no dispersed
phase in solid) both exactly 0; $\Delta p$ relative error
$2.60\times10^{-3}$, $\max|\mathbf{u}| = 1.22\times10^{-4}$; 592 leaves
over levels 4–7. The negative control fails on P1 at $i = 0$ with a
$6.25\times10^{-4}$ deviation in `fm`.

`staticFilmTube` — $N = 256$, $\Delta = 7.81\times10^{-3}$,
$t_{end} = 2$, production ratios $\rho_d/\rho_c = 10^{-3}$,
$\mu_d/\mu_c = 10^{-2}$. Sustained $\max|\mathbf{u}|$ is the peak over
$t \ge t_{end}/2$:

| gap [cells] | `TOLERANCE` | rel. $\Delta p$ error | sustained $\max\|\mathbf{u}\|$ |
|---|---|---|---|
| 4 | $10^{-5}$ | $5.36\times10^{-5}$ | $1.53\times10^{-4}$ |
| 6 | $10^{-5}$ | $8.36\times10^{-5}$ | $1.66\times10^{-4}$ |
| 8 | $10^{-5}$ | $1.01\times10^{-4}$ | $1.59\times10^{-4}$ |
| 4 | $10^{-6}$ | $5.22\times10^{-5}$ | $8.19\times10^{-6}$ |
| 4 | $10^{-4}$ (production) | $\sim9\times10^{-5}$ | $1.2$–$1.8\times10^{-3}$ |

## What the two sweeps measure

**Pressure-jump error does not degrade as the gap narrows.** At $N=256$
the relative $\Delta p$ error at a 4-cell gap is $5.4\times10^{-5}$, some
370× inside the 2% gate, and it is slightly *larger* at 8 cells. Across
grids it falls by ~6× from $N=128$ to $N=256$. The error tracks
$R_d/\Delta$, not wall proximity, so the non-embed-aware height-function
stencil shows no detectable damage down to 4 cells of separation at
production property ratios. The ordering of the three gaps is *not*
resolved at $N=256$: the differences are comparable to the within-run
temporal fluctuation of $\Delta p$.

**The spurious-current amplitude is set by the projection tolerance, not
by curvature, once the density ratio is large.** In `staticFilmTube`,
tightening `TOLERANCE` from $10^{-5}$ to $10^{-6}$ moves $\Delta p$ by
$1.3\times10^{-6}$ relative — immaterial, so the pressure-jump error is
genuine discretisation error — but moves $\max|\mathbf{u}|$ by a factor
of about 20. The tabulated velocities are therefore an upper bound on the
curvature-driven current, not a measurement of it. Record the
`TOLERANCE` any reported velocity was measured at.

This does **not** transfer to the matched-property cases. Repeating the
ladder on `laplaceEmbedTube` and `laplaceEmbedTubeAdapt` at $10^{-6}$
moves their $\max|\mathbf{u}|$ by 0.1–1.3%, so those amplitudes are
curvature-driven. Note that N-convergence alone cannot distinguish the
two: the demonstrably tolerance-set amplitude in `staticFilmTube` also
falls by ~3.5× from $N=128$ to $N=256$. Only the ladder separates them.

**Reproducibility.** Three independent 16-thread runs of
`staticFilmTube` spread by 0.1% at $10^{-6}$ and 9% at $10^{-5}$, but by
45% at $10^{-4}$ ($1.25$, $1.74$ and $1.81\times10^{-3}$) — itself
consistent with a residual-dominated amplitude. The $10^{-5}$ and
$10^{-6}$ rows above are one run; treat the $10^{-4}$ row as an order of
magnitude only.

## Basilisk version and the axisymmetric embed advection patch

Results here were produced with Basilisk at patch **2026-07-03**.

On 2026-08-03 two upstream patches fixed a defect in `src/embed.h`
affecting `embed.h` combined with `axi.h`
([changes](https://basilisk.fr/src/?changes=20260803152337),
[test](https://basilisk.fr/src/test/pipe-axi-embed.c)). In
`update_tracer()` the flux divergence was divided by `Delta` rather than
`Delta*cm[]`; under the axisymmetric metric `cm = y*cs`, so the
advection term acquired a spurious factor `y`. The reported symptom is a
laminar pipe entrance length roughly six times too short. The developed
Poiseuille state is unaffected.

That combination is exactly the one this solver uses, so it was checked
rather than assumed. `update_tracer()` is reached only from
`advection()` in `bcg.h`, which `navier-stokes/centered.h` calls only
inside `if (!stokes)`. `navier-stokes/conserving.h` sets `stokes = true`
in its `defaults` event, so with momentum-conserving VOF advection that
branch is never taken. Confirmed at runtime: a breakpoint on
`update_tracer` was never hit over a complete run of
`simulationCases/bretherton.c`, which then exited normally. The function
is compiled into the binary but not executed.

Independently, these runs are viscous-dominated: with `La = 1` the
Reynolds number is `La*Ca`, so 0.005 to 0.07 across the campaign, and the
corrupted term is the inertial one.

Both points should be rechecked if the solver ever drops
`conserving.h`, or runs at a Reynolds number where inertia matters.

## What this suite does not establish

- Nothing dynamic is tested. Every case has a quiescent or
  prescribed-advection exact solution. The Bretherton film is non-spherical, carries a
  lubrication pressure gradient and is advected; none of that is tested.
- Only one case adapts, and it runs matched properties with the drop far
  from the wall. `staticFilmTube` puts the interface near the wall but on
  a uniform grid. Nothing tests the two together, which is the
  combination production actually runs.
- The negative control is weak. It fires on a single face at the
  domain boundary where the flux is zero regardless, and P3 is
  bit-identical between the guarded and guard-less builds. It shows the
  probe detects metric staleness; it does not show that omitting the sync
  would corrupt a production run. P2 never fires, because no case puts
  the interface close enough to the wall on an adapted grid.
- No case except `staticFilmTube` carries a property contrast. With
  $\rho_1 = \rho_2$ and $\mu_1 = \mu_2$ the property construction in
  `two-phase-generic.h` is constant, so the matched-property cases cannot
  detect a regression in that path at all.
- `staticFilmTube`'s 8-cell gap had not plateaued by $t_{end}$:
  the last-quarter/second-quarter mean ratio was 1.51 and 1.63 in two
  independent runs, both above the 1.05 threshold the case gates on. Its
  $\max|\mathbf{u}|$ is therefore a lower bound, and the case prints
  that advisory.
