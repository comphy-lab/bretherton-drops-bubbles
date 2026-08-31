# Bretherton drops and bubbles

Axisymmetric Basilisk solver for confined drops and bubbles in capillary
tubes: embedded tube wall (`embed.h`), VOF interface (`two-phase.h`),
surface tension, Bretherton validation.

## Structure

```
bretherton-drops-bubbles/
├── src-local/         # project headers: parameter accessors + embed/VOF compatibility
├── simulationCases/   # Basilisk entry points; numbered case output directories (gitignored)
├── testCases/         # smoke + verification cases and their driver
├── postProcess/       # snapshot and log analysis
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
- Smoke-test before committing solver changes:
  `bash testCases/runSmokeTests.sh`. Tests are classified by evidence
  source (see `testCases/README.md`): smoke = compiles and runs a few
  steps; verification = exact solutions; validation = Bretherton /
  Taylor data via the root sweep files on production hardware.
- Never hardcode a machine-local `qcc` path; resolve via `PATH` or
  `.project_config` (gitignored).
- Do not commit `basilisk/`, `.comphy-basilisk`, case outputs, or
  `CLAUDE.md` (gitignored; contains only `@AGENTS.md`).
- Component READMEs and `docs/` are public-candidate. Record
  simulation-time numbers internally first; promotion of findings into
  the README requires explicit approval.
- Use `publication-plots` for every figure, including diagnostics.
