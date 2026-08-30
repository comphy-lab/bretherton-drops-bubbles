# Bretherton drops and bubbles

Axisymmetric Basilisk solver for long drops and bubbles in capillary tubes, aimed at recovering the classical Bretherton film-thickness law and then comparing drop and bubble counterparts.

## Overview

A long gas bubble or liquid drop moving through a liquid-filled tube leaves a thin lubricating film on the wall. For a bubble at small capillary number, lubrication theory gives a film thickness that scales as \(\mathrm{Ca}^{2/3}\). This repository will hold the DNS that recover that law and then extend it to density- and viscosity-ratio effects for drops.

The solver has not been scaffolded yet. Cases will later live in `simulationCases/` with shared headers in `src-local/` and analysis in `postProcess/`.

## Requirements

- Basilisk C (`qcc`)
- MPI for production runs, when the solver is added

## License

See the repository license once one is added.
