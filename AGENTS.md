# Bretherton drops and bubbles

Axisymmetric Basilisk solver for confined drops and bubbles in tubes. Recover the classical Bretherton film-thickness law before adding drop/bubble extensions.

## Structure

The repository is currently an empty scaffold. When the solver is added, use the CoMPhy Basilisk layout:

```
bretherton-drops-bubbles/
├── src-local/         # project headers and parameter accessors
├── simulationCases/   # Basilisk entry points and numbered case directories
├── postProcess/       # analysis scripts
├── default.params     # runtime defaults
└── sweep.params       # parameter-sweep contract
```

## Building

```bash
# Compile with qcc once a case exists
qcc -O2 -Wall simulationCases/<case>.c -o simulation -lm
```

## Guidelines

- Units: dimensionless, with capillary number \(\mathrm{Ca}\) as the primary control parameter.
- Do not invent solver files; wait for an explicit scaffold request.
- Component READMEs and `docs/` are public-candidate. Record simulation-time numbers internally first.
- Use `publication-plots` for every figure, including diagnostics.
- Local `CLAUDE.md` is gitignored and must contain only `@AGENTS.md`.
