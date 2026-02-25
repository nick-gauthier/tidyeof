# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

tidyeof is an R package for empirical orthogonal function (EOF) analysis in the tidyverse framework. It extracts modes of variability from spatiotemporal data (stars objects), runs diagnostics, and performs spatial downscaling via canonical correlation analysis (CCA).

## Common Commands

```r
devtools::load_all()       # Load package for interactive development
devtools::test()           # Run full test suite
devtools::test(filter = "get_patterns")  # Run a single test file
devtools::check()          # Full R CMD check
devtools::document()       # Regenerate man/ pages and NAMESPACE from roxygen2
```

Vignettes use Quarto (.qmd), not R Markdown.

## Architecture

### S3 Class Hierarchy

Three main S3 classes, each with print/plot/predict/summary/subset methods:

- **`patterns`** — output of `patterns()`. Contains eigenvalues, EOFs (spatial loadings), and amplitudes (time series). Created in `patterns_core.R`, class defined in `patterns_class.R`.
- **`coupled_patterns`** — output of `couple()`. CCA results pairing predictor and predictand patterns. Defined in `couple.R`.
- **`common_patterns`** — output of `common_patterns()`. Joint EOF across multiple datasets with per-source amplitudes. Defined in `common_patterns.R`.

### Core Pipeline

1. **Climatology** (`get_climatology.R`) — compute monthly climatology and anomalies from stars objects
2. **EOF extraction** (`patterns_core.R`) — PCA on the space-time matrix; auto-selects `irlba::prcomp_irlba()` for large datasets (>500k elements) vs base `prcomp()`
3. **Coupling** (`couple.R`) — CCA-based downscaling linking coarse predictor EOFs to fine predictand EOFs
4. **Cross-validation** (`tune_cca.R`, `prep_cv.R`) — contiguous temporal folds for tuning EOF/CCA hyperparameters

### Spatial Data Handling

The package works with `stars` objects and supports two spatial formats:
- **Raster grids** (x, y, time dimensions)
- **sf geometry** (geometry, time dimensions)

Auto-detection via `has_geometry_dimension()` and `get_spatial_dimensions()` in `stars_utils.R`. CRS and units are preserved throughout.

## Code Conventions

- Tidyverse style; tidy evaluation throughout
- User-facing messages via `cli`; errors via `rlang::abort()` with custom `tidyeof_*` error classes
- Documentation via roxygen2 (`@export`, `@keywords internal`)
- Test data: single file `inst/testdata/prism_test.RDS` (51×51 grid, 36 monthly timesteps of PRISM temperature in °C), loaded globally in `tests/testthat.R` as `prism`
