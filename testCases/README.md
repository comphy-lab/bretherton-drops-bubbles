# Test cases

Basilisk tests are not unit tests: every case here proves, at minimum,
that the code compiles and integrates a few time steps (a *smoke test*).
Cases are classified by their evidence source:

| Case | Class | Evidence |
|------|-------|----------|
| `embedAxiVofAdvection.c` | verification | exact solution $r(t)=\sqrt{r_0^2+2t}$ of axisymmetric VOF advection with an embedded wall (adapted from upstream `src/test/missing_metric.c`) |
| `laplaceEmbedTube.c` | verification | Young–Laplace pressure jump $\Delta p = 2\sigma/R_d$ and zero velocity for a static droplet inside an embedded tube |
| `smoke.params` + `simulationCases/bretherton.c` | smoke test | compiles and integrates ~60 steps at coarse resolution; checks the log is produced and nothing blows up |

The *validation* campaign (comparison against Bretherton's law and
Taylor/Aussillous–Quéré data, i.e. evidence independent of the
implemented equations) is configured by `sweep.params` (bubble) and
`sweep-drop.params` (drop) at the repository root and is meant for
production hardware, not this directory.

## Running

```bash
bash testCases/runSmokeTests.sh              # everything
bash testCases/runSmokeTests.sh --skip-main  # only the two verification cases
```

Both verification cases print `PASS`/`FAIL` and their measured errors.
Reference results on a single core (Basilisk `qcc`, quadtree):

- `embedAxiVofAdvection`: max relative interface-position error
  ~3×10⁻⁵ (tolerance 10⁻²).
- `laplaceEmbedTube`: pressure-jump error ~6×10⁻⁴ (tolerance 2×10⁻²),
  spurious currents ~4×10⁻⁵ σ/μ (tolerance 10⁻³) at t = 1.

Build artefacts go to `testCases/build-*/` and the smoke case to
`simulationCases/9999/`; both are gitignored.
