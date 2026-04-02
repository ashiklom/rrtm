## Why

The repository documents a Bayesian inversion workflow for PROSPECT in `vignettes/inversion.Rmd`, but there is no parallel prototype showing how to express that model in Nimble. Evaluating Nimble is easier if the repository includes both a minimal integration attempt that calls `rrtm::prospectd` directly and a deeper prototype that rewrites the forward model in `nimbleCode`.

## What Changes

- Add prototype implementations under `misc/` that reproduce the inversion workflow described in `vignettes/inversion.Rmd` using Nimble in two ways: a direct-call variant and a `nimbleCode` rewrite variant.
- Ensure the prototype material is kept out of the package build by adding `misc/` to `.Rbuildignore`.
- Document the inputs, parameterization, execution steps, and trade-offs of the two prototype paths so they can be compared for future development.

## Capabilities

### New Capabilities
- `nimble-inversion-prototype`: Provide non-package prototypes that run the inversion workflow from `vignettes/inversion.Rmd` with Nimble, covering both a direct-call integration and a `nimbleCode` rewrite while keeping the work isolated from package builds.

### Modified Capabilities

None.

## Impact

- Affects repository-level build exclusions via `.Rbuildignore`.
- Adds prototype R code and supporting notes under `misc/`.
- Introduces two Nimble-based workflows: one that depends directly on `rrtm::prospectd` behavior and one that re-expresses that behavior in `nimbleCode`, without changing package APIs.
