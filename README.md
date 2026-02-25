# tidyeof

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

An R package for conducting empirical orthogonal function (EOF) analysis in the tidyverse framework. Functions to isolate modes of variability from spatiotemporal data, run various diagnostics, and use these patterns for spatial downscaling.

## Installation

```r
# Install from source
devtools::install_github("nick-gauthier/tidyeof")
```

## Usage

```r
library(tidyeof)

# Extract EOF patterns from spatiotemporal data
pat <- patterns(climate_data, k = 5)

# View patterns
plot(pat)
screeplot(pat)

# Project new data onto existing patterns
amplitudes <- project_patterns(pat, new_data)

# Reconstruct fields from patterns
reconstructed <- reconstruct(pat, amplitudes)
```

## Pattern Coupling for Downscaling

```r
# Extract patterns from predictor and response fields
pred_patterns <- patterns(coarse_data, k = 5)
resp_patterns <- patterns(fine_data, k = 5)

# Couple patterns using CCA
coupled <- couple(pred_patterns, resp_patterns, k = 3)

# Make predictions on new data
downscaled <- predict(coupled, new_coarse_data)
```

## Cross-Validation

```r
# Prepare CV folds (expensive, run once)
cv <- prep_cv_folds(coarse_data, fine_data,
                    kfolds = 5, max_k_pred = 10, max_k_resp = 10)

# Tune hyperparameters
results <- tune_cca(cv, k_pred = 1:10, k_resp = 1:10)

# Summarize results
summary <- summarize_cv(results, metric = "rmse", minimize = TRUE)
```
