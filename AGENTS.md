# Bretherton drops and bubbles

Axisymmetric Basilisk solver for confined drops and bubbles in capillary
tubes: embedded tube wall (`embed.h`), VOF interface (`two-phase.h`),
surface tension, Bretherton validation.

## Structure

```
bretherton-drops-bubbles/
├── src-local/         # project headers: parameter accessors + embed/VOF compatibility
├── simulationCases/   # Basilisk entry points; numbered case output directories (gitignored)
├── verificationCases/ # exact-solution cases and their driver
├── testCases/         # smoke test only
├── postProcess/       # snapshot and log analysis
├── runTests.sh        # entry point: verification cases, then smoke test
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
  separation is deliberate: `testCases/` = smoke, compiles and runs a few
  steps against no comparator; `verificationCases/` = exact solutions of
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
  `KAPPA` criteria plus the explicit `refine()` in the init event. Do not
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
