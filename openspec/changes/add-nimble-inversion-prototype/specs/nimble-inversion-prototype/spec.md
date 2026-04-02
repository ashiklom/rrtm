## ADDED Requirements

### Requirement: Provide standalone Nimble inversion prototypes
The repository SHALL include standalone prototypes under `misc/` that reproduce the Bayesian inversion workflow described in `vignettes/inversion.Rmd` using Nimble while leaving the installed package API unchanged.

#### Scenario: Prototype location is isolated from package code
- **WHEN** a contributor inspects the implementation for the Nimble inversion workflow
- **THEN** the prototype files are located under `misc/` rather than `R/`, `src/`, or `vignettes/`

#### Scenario: Both implementation strategies are present
- **WHEN** a contributor inspects the Nimble prototype work
- **THEN** the repository provides one variant that attempts to call `rrtm::prospectd` directly and one variant that rewrites the forward model in `nimbleCode`

#### Scenario: Prototype follows the documented inversion flow
- **WHEN** a contributor reads the prototype implementation
- **THEN** each variant includes the same major stages as the vignette: loading spectral observations, interpolating to PROSPECT wavelengths, specifying priors and likelihood around the PROSPECT forward model, and running a Nimble-based Bayesian fit

### Requirement: Preserve the vignette model assumptions in both prototype variants
The Nimble prototypes SHALL preserve the inversion parameter set and observation model from `vignettes/inversion.Rmd`, including inference for `N`, `Cab`, `Car`, `Cbrown`, `Canth`, `Cw`, `Cm`, and residual standard deviation.

#### Scenario: Parameter set matches the vignette
- **WHEN** the prototype defines its latent parameters
- **THEN** it includes the seven PROSPECT parameters plus an inferred residual scale parameter

#### Scenario: Invalid forward-model states are handled explicitly
- **WHEN** either prototype variant fails or produces illegal reflectance values during evaluation
- **THEN** that variant applies a defined guard path that prevents those states from being treated as valid high-probability samples

#### Scenario: Rewrite variant remains comparable to the package model
- **WHEN** the `nimbleCode` rewrite variant is exercised on the inversion workflow
- **THEN** the implementation includes a documented way to compare its forward-model outputs against `rrtm::prospectd` for shared parameter inputs

### Requirement: Exclude prototype artifacts from package builds
The repository SHALL ignore the `misc/` directory during package builds so experimental Nimble prototype files are not bundled as package build artifacts.

#### Scenario: Build ignore includes prototype directory
- **WHEN** `.Rbuildignore` is evaluated during package build preparation
- **THEN** entries in `misc/` are excluded from the build input set
