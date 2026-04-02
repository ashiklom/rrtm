## 1. Repository Setup

- [x] 1.1 Create the `misc/` directory structure for the Nimble prototype and add `^misc/` to `.Rbuildignore`.
- [x] 1.2 Add a short prototype-facing README or script header that explains required optional packages, the two implementation variants, and how to run them outside the package build.

## 2. Shared Prototype Foundations

- [x] 2.1 Port the observation loading and interpolation steps from `vignettes/inversion.Rmd` into a standalone `misc/` script.
- [x] 2.2 Define shared priors, parameter naming, and observation-model conventions so both Nimble variants match the vignette as closely as possible.

## 3. Direct-Call Variant

- [x] 3.1 Implement the Nimble-compatible forward-model wrapper around `rrtm::prospectd`, including explicit handling for runtime errors and illegal reflectance values.
- [x] 3.2 Define and run the direct-call Nimble model using the vignette parameter set and posterior predictive summaries.
- [x] 3.3 Record any limitations, workarounds, or incompatibilities encountered when trying to call `rrtm::prospectd` directly from Nimble.

## 4. `nimbleCode` Rewrite Variant

- [x] 4.1 Re-express the required `rrtm::prospectd` forward-model calculations in `nimbleCode` for the inversion workflow.
- [x] 4.2 Define and run the rewrite-based Nimble model using the same parameter set, priors, and posterior predictive summaries.
- [x] 4.3 Add a comparison step that checks rewrite outputs against `rrtm::prospectd` for shared parameter inputs and documents discrepancies.

## 5. Prototype Verification

- [x] 5.1 Run both prototype variants far enough to verify that their Nimble models initialize without immediate model-definition errors.
- [x] 5.2 Confirm that the prototype remains outside package build artifacts by checking the `misc/` exclusion in build-related files.
- [x] 5.3 Summarize the relative strengths, limitations, and follow-up work for the direct-call and rewrite variants.
