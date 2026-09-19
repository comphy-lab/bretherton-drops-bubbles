# Bretherton drops and bubbles

Axisymmetric Basilisk solver for confined drops and bubbles in capillary
tubes: embedded tube wall (`embed.h`), VOF interface (`two-phase.h`),
surface tension, Bretherton validation.

## Structure

```
bretherton-drops-bubbles/
├── src-local/         # project headers: parameter accessors, embed/VOF compatibility,
│                      # central-film observation contract
├── simulationCases/   # Basilisk entry points (lab frame bretherton.c, co-moving
│                      # bretherton-comoving.c); numbered case output directories (gitignored)
├── verificationCases/ # exact-solution cases and their driver
├── testCases/         # observer unit tests and execution/restart checks
├── postProcess/       # snapshot and log analysis (`central_film.py` replays observation windows from logs)
├── runContinuation.py # Ca continuation ladder on top of runSimulation.sh (co-moving case)
├── continuation/      # station parameter files for the co-moving ladder
├── runTests.sh        # entry point: unit tests, verification cases, then smoke tests
├── default.params     # runtime defaults (bubble)
├── sweep.params       # bubble validation sweep contract
└── sweep-drop.params  # drop counterpart sweep contract
```

## Non-negotiables

- Units: dimensionless with repeating variables $\sigma$, $\mu_c$ and the
  volume-equivalent radius $R$. Velocities are capillary numbers; do not
  reintroduce dimensional parameters.
- Every `solid()` call and every `adapt_wavelet()` must be followed by
  `embed_axi_metric_sync()` and (after adaptation) `vof_solid_cleanup(f)`
  from `src-local/embed-vof-tube.h`. Dropping either silently breaks the
  axi+embed metric or the near-wall curvature.
- The interface must never touch the embedded wall: keep the film
  resolved by at least ~4 cells at `MAXlevel` and watch the film warning
  in the case log.
- The time-step cap is `DT` (runtime key `dtmax` maps onto it);
  assigning the Basilisk global `dtmax` directly is a no-op because
  `centered.h` resets it every iteration.
- Event-loop increments (`t += tsnap`) are classified before `main()`
  runs: keep static non-zero initialisers on `tsnap`-like globals.
- One parameter pathway only: `src-local/parse_params.h` + `params.h`
  read `case.params` (argv[1]); the shell runners use
  `get_param_value`/`set_param_in_file`. Do not add parallel parsers.
- Run cases through `runSimulation.sh` / `runParameterSweep.sh`; they
  create `simulationCases/<CaseNo>/` (CaseNo >= 1000; 9999 is reserved
  for the smoke test).
- Run the evidence suite before committing solver changes:
  `bash runTests.sh`. Cases are separated by evidence source, and the
  separation is deliberate: `testCases/` = software unit, execution and
  transfer checks; `verificationCases/` = exact solutions of
  the implemented equations (see `verificationCases/README.md`);
  validation = independent experimental data (Taylor 1961, Aussillous &
  Quéré 2000) via the root sweep files on production hardware. Bretherton
  (1961) is an asymptotic solution of the lubrication limit, not
  independent data — comparing against it is a limiting-case check, not
  validation.
- `adapt_wavelet()` on `cs` does **not** refine a flat embedded wall.
  `embed.h` sets `cs.prolongation = fraction_refine`, which is exact for
  any planar interface, so the tube wall carries identically zero wavelet
  error at any `csErr`. Wall and film refinement come from the `f` and
  `KAPPA` criteria plus explicit initialization and measured-film refinement.
  `filmMinLevel` floors the interface and wall-facing liquid band;
  `filmCells` uses a measured thickness and is capped by `MAXlevel`. Do not
  add a `csErr` and assume the wall is resolved; check the cut-cell level
  range instead.
- A reported velocity is meaningless without the `TOLERANCE` it was
  measured at. `verificationCases/staticFilmTube.c` shows that once the
  density ratio is large, the spurious-current amplitude in a static
  configuration is set by the projection residual rather than by
  curvature: a 10x tighter tolerance moves it ~19x while leaving the
  pressure jump unchanged. Always record the tolerance alongside any
  velocity, and run the ladder before attributing a current to the
  surface-tension discretisation.
- Never hardcode a machine-local `qcc` path; resolve via `PATH` or
  `.project_config` (gitignored).
- Do not commit `basilisk/`, `.comphy-basilisk`, case outputs, or
  `CLAUDE.md` (gitignored; contains only `@AGENTS.md`).
- Component READMEs and `docs/` are public-candidate. Record
  simulation-time numbers internally first; promotion of findings into
  the README requires explicit approval.
- Use `publication-plots` for every figure, including diagnostics.

## Co-moving window (`simulationCases/bretherton-comoving.c`)

- The frame translates at speed `U` set by a critically damped controller on
  the VOF centroid (`src-local/comoving-frame.h`). Every velocity boundary
  value is shifted by `-U` and the frame acceleration enters as the uniform
  body force `-dU/dt` through the acceleration face vector. Do not add a
  density weight: Basilisk's `a` is an acceleration per unit mass and the
  transformation is exact for both phases.
- On the embedded wall `u.n[embed]` is the axial component and `u.t[embed]`
  the radial one (Basilisk convention, `src/test/couette.c`). The wall moves
  at `-U`; its geometry never moves.
- In the frame the mean liquid speed is `Ca_in - U_b < 0`: liquid enters
  through the front end and leaves through the rear. With Poiseuille minus `U`
  imposed at both ends the pressure floats; report gauges against the mean
  front-end pressure. Both end columns are held at `padLevel` so the two
  discrete boundary fluxes match; the residual is logged as `Qrear`/`Qfront`.
- Steady state is judged in film renewals `L_b/U`, not in advance windows.
  The film freshness gate is `freshFront = xTipF - displacement`.
- Frame scalars are not in a dump. A same-case resume reads `frame-state`; a
  continuation station receives `Uframe0`, `xTarget` and `CaPrev` through
  `case.params` (`runContinuation.py` writes them from the seed's
  `station.json`) and must start with exactly the seed frame speed, because
  the dumped velocity field lives in that frame.
- One `Ca` per case directory is now the continuation contract as well as
  the log-schema rule; never resume a directory with a different `Ca`.

