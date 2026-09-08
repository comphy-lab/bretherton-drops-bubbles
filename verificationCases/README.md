# Verification cases

These programs compare the implemented model with analytical solutions or
exact geometric integrals. Experimental film-thickness comparisons belong in
[the Newtonian validation report](../docs/Newtonian-Validation/).

| Case | Check |
|---|---|
| `embedAxiVofAdvection.c` | Prescribed axisymmetric interface advection |
| `laplaceEmbedTube.c` | Young–Laplace pressure jump on a sequence of uniform grids |
| `laplaceEmbedTubeAdapt.c` | Young–Laplace pressure jump and coupling checks on an adaptive grid |
| `staticFilmTube.c` | Static interface near the wall, with density and viscosity contrasts |
| `tubeGeometry.c` | Exact tube volume, axial area and constant-flow flux residual through repeated adaptation |

Use the project-local Basilisk installation configured by `.project_config`.
From the repository root:

```bash
bash runTests.sh
bash runTests.sh --verification
bash verificationCases/runVerification.sh verificationCases/tubeGeometry.c
VERIFICATION_THREADS=16 bash runTests.sh
```

Builds and diagnostic output are written to ignored `build-<case>/`
directories. Each verification program returns nonzero on failure and prints
a bare `PASS` or `FAIL` line. Record the toolchain and thread count when
retaining results.

The driver also builds `laplaceEmbedTubeAdapt.c` with
`-DSKIP_EMBED_GUARDS`. This negative control uses active-cell metric
initialisation without the complete stored-tree tube reconstruction, then
omits the post-adaptation sync and cleanup calls. It must fail a consistency
probe. It tests the combined compatibility pathway rather than isolating
the contribution of each call.
