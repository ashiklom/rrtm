# Nimble Inversion Prototype

This directory contains a standalone prototype of the Bayesian inversion workflow from `vignettes/inversion.Rmd`.
It is intentionally kept outside the package build and is excluded by `.Rbuildignore`.

## Requirements

- `nimble`
- `rrtm`
- `pkgload` if you want the script to load the checkout in place instead of an installed `rrtm`

## Variants

- `direct`: uses `nimbleRcall` to evaluate `rrtm::prospectd` from inside the Nimble model
- `rewrite`: re-expresses the `prospectd()` and `gpm()` calculations in native `nimbleFunction` code, using the package lookup tables as explicit constants

## Run

From the repository root:

```r
Rscript misc/nimble-inversion/prototype.R --variant both --iterations 200 --burnin 100
```

Useful options:

- `--variant direct`
- `--variant rewrite`
- `--summary-samples 25`

The script loads the NASA FFT example used by the vignette, interpolates it to the PROSPECT wavelengths, runs the requested Nimble variant, and prints compact summaries.

## Notes

- The direct-call path is the smallest change from the current vignette, but every likelihood evaluation crosses back into ordinary R through `nimbleRcall`, so it is mainly a feasibility test.
- The rewrite path stays inside Nimble for the forward model, but it is substantially more code and can diverge numerically from `rrtm::prospectd` if any part of the port is imperfect.
- The script compares rewrite outputs against `rrtm::prospectd` on shared parameter inputs and reports the discrepancy so the two approaches can be judged side by side.

## Current Findings

- Both variants now build, compile, and run short MCMC chains from the prototype script.
- The rewrite variant currently matches the direct/package forward model exactly for the tested parameter comparisons reported by the script (`max_abs_diff = 0`, `rmse = 0`).
- The direct variant reuses `rrtm::prospectd`, so it inherits warnings from the package forward model under some high-absorption proposals; the prototype now muffles those repeated warnings to keep MCMC output readable.
- The rewrite variant required a small Nimble-specific coding adjustment: scalar `||` conditions in the native rewrite generated invalid C++, so the script uses separate scalar `if` checks instead.
- These runs are still prototype-scale checks, not convergence-validated inference runs.
