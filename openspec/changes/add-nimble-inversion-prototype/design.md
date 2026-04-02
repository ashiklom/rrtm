## Context

`vignettes/inversion.Rmd` currently demonstrates two inversion approaches around `rrtm::prospectd`: a direct least-squares fit via `optim()` and a Bayesian fit via `BayesianTools`. The requested change is a prototype, not a package feature, so the implementation should live outside `R/`, `man/`, and the package build path while still being easy to run and inspect. The repository does not currently contain a `misc/` directory, and `DESCRIPTION` does not list `nimble` as a package dependency, so the prototype should be documented as an optional, analyst-facing workflow rather than integrated package functionality. The revised goal is to compare two Nimble approaches: one that tries to call `rrtm::prospectd` directly from Nimble machinery and one that rewrites the forward model in `nimbleCode`.

## Goals / Non-Goals

**Goals:**
- Re-express the Bayesian inversion example from `vignettes/inversion.Rmd` as a Nimble-based prototype.
- Compare two implementation strategies for that prototype: direct invocation of `rrtm::prospectd` and a native `nimbleCode` rewrite of the forward model.
- Keep the prototype isolated in `misc/` and excluded from `R CMD build`.
- Preserve the core modeling assumptions from the vignette: the same observed spectrum preprocessing, the same PROSPECT-D forward model, similar priors, and a Gaussian observation model with residual standard deviation as an inferred parameter.
- Make the prototype understandable enough to serve as a starting point for future package-grade work.

**Non-Goals:**
- Adding Nimble as a required package dependency for `rrtm`.
- Refactoring `rrtm::prospectd` or exposing new public APIs.
- Guaranteeing production-quality runtime performance, convergence diagnostics, or package-level tests for the prototype.

## Decisions

### Build the prototype as a standalone R script under `misc/`
The implementation will be a standalone script, with any minimal supporting notes colocated in `misc/`, because the request is for a prototype and the repository currently treats experimental material as non-package content. This avoids accidental API commitments and keeps the package namespace unchanged.

Alternative considered: adding the code as a vignette or exported helper. Rejected because that would imply stronger support guarantees, create a package dependency question for Nimble, and make the experimental workflow part of the published package surface.

### Implement two prototype variants under `misc/`
The change should produce two prototype variants: one variant that attempts to integrate Nimble with the existing `rrtm::prospectd` implementation, and a second variant that rewrites the forward model in `nimbleCode`. Keeping both variants side by side allows the repository to compare integration feasibility, code complexity, and likely performance trade-offs.

Alternative considered: choosing only one approach up front. Rejected because the user explicitly wants both paths explored and the comparison is part of the value of the prototype.

### Keep one variant as a thin wrapper around `rrtm::prospectd`
The first variant should call `rrtm::prospectd` from a custom Nimble-compatible function rather than rewriting PROSPECT-D. This provides the lowest-effort path to test whether Nimble can work with the existing forward model closely enough to support inference.

Alternative considered: dropping the direct-call path once the rewrite path exists. Rejected because the direct-call attempt is useful even if it reveals limitations; those limitations are part of the intended outcome.

### Build a second variant as a `nimbleCode` rewrite of the forward model
The second variant should re-express the forward model logic in `nimbleCode` so the entire likelihood path is native to Nimble. This provides a clearer picture of what a deeper Nimble integration would require and may avoid some constraints of calling ordinary R code during sampling.

Alternative considered: stopping at partial pseudocode for the rewrite. Rejected because the prototype is meant to compare implementation strategies in code, not only on paper.

### Mirror the vignette priors and likelihood as closely as Nimble allows
The Nimble model should infer `N`, `Cab`, `Car`, `Cbrown`, `Canth`, `Cw`, `Cm`, and `rsd`, using the same distribution families and parameter interpretation described in the vignette. If exact parameterization differs between packages, the prototype should document the mapping explicitly in comments.

Alternative considered: simplifying the priors to generic weakly informative normals. Rejected because the comparison to the existing vignette would become less meaningful.

### Handle invalid forward-model evaluations with explicit safeguards in both variants
Because the vignette already notes that some parameter combinations can produce invalid spectra or runtime errors, each prototype variant should include a deterministic guard path that returns a poor likelihood contribution when the forward model fails or yields illegal reflectance values.

Alternative considered: letting sampler failures propagate. Rejected because that would make the prototype brittle and hard to evaluate.

### Exclude `misc/` from package builds via `.Rbuildignore`
The repository already excludes `data-raw/`, `extdata/`, and `scratch/`; `misc/` should be handled similarly so prototype code does not appear in package builds or CRAN-facing artifacts.

Alternative considered: leaving `misc/` tracked but build-visible. Rejected because it conflicts with the user requirement and blurs the line between prototype and package content.

## Risks / Trade-offs

- [Nimble integration with `rrtm::prospectd` may require custom function registration details or may not work cleanly] -> Mitigation: keep the direct-call path as an explicit experiment and document the exact limitations encountered.
- [A `nimbleCode` rewrite of PROSPECT-D may be substantial and can diverge numerically from `rrtm::prospectd`] -> Mitigation: structure the rewrite to compare outputs against `rrtm::prospectd` on shared parameter inputs and record any approximation or omission.
- [Running PROSPECT inside MCMC can be slow] -> Mitigation: frame the implementation as a prototype, start with a modest iteration count, and keep data handling simple.
- [Nimble distribution parameterization may differ from `distributions3`] -> Mitigation: annotate the chosen parameter mapping in the script and keep priors aligned by family and scale.
- [Prototype code can drift from the vignette workflow over time] -> Mitigation: structure both variants in the same stages as the vignette: load data, interpolate observations, define model, run MCMC, summarize predictions.
