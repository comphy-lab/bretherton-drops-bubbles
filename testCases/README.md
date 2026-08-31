# Smoke test

This directory holds the repository's smoke test and nothing else.

A smoke test proves that `simulationCases/bretherton.c` compiles and
integrates a few dozen time steps at coarse resolution without producing
a truncated log or blowing up. It compares against **nothing**, so it
supports no claim about the discretisation or about the physics. Its
value is that it fails loudly and in seconds when a change breaks the
build or the parameter path.

| File | Role |
|------|------|
| `smoke.params` | coarse, short-run parameters (`MAXlevel = 9`, `tmax = 0.05`) |
| `runSmokeTests.sh` | builds and runs the case, then checks the log |

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
bash runTests.sh                  # verification cases, then this
```

The build goes to `simulationCases/9999/`, which is gitignored.
