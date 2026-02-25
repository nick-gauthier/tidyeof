# Package index

## Core EOF Analysis

Extract and work with EOF patterns

- [`patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/patterns.md)
  : Get EOFs and PCs from spatiotemporal data
- [`project_patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/project_patterns.md)
  : Project New Data onto Existing Patterns
- [`project_patterns_multiple()`](https://nick-gauthier.github.io/tidyEOF/reference/project_patterns_multiple.md)
  : Project Multiple Datasets onto Patterns
- [`reconstruct()`](https://nick-gauthier.github.io/tidyEOF/reference/reconstruct.md)
  : Reconstruct spatial field from EOF amplitudes
- [`area_weights()`](https://nick-gauthier.github.io/tidyEOF/reference/area_weights.md)
  : Compute area-based weights for spatial data

## Climatology

Compute and restore climatological baselines

- [`get_climatology()`](https://nick-gauthier.github.io/tidyEOF/reference/get_climatology.md)
  : Calculate climatological mean and standard deviation for spatial
  data
- [`get_anomalies()`](https://nick-gauthier.github.io/tidyEOF/reference/get_anomalies.md)
  : Calculate anomalies from a climatological mean
- [`restore_climatology()`](https://nick-gauthier.github.io/tidyEOF/reference/restore_climatology.md)
  : Restore original field from anomalies and climatology

## CCA Downscaling

Couple patterns and make predictions

- [`couple()`](https://nick-gauthier.github.io/tidyEOF/reference/couple.md)
  : Couple Pattern Relationships Using CCA
- [`get_canonical_correlations()`](https://nick-gauthier.github.io/tidyEOF/reference/get_canonical_correlations.md)
  : Get Canonical Correlations from Coupled Patterns
- [`get_canonical_patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/get_canonical_patterns.md)
  : Get Canonical Spatial Patterns from Coupled Patterns
- [`get_canonical_variables()`](https://nick-gauthier.github.io/tidyEOF/reference/get_canonical_variables.md)
  : Get Canonical Variables from Coupled Patterns

## Common Patterns

Joint EOF analysis across datasets

- [`common_patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/common_patterns.md)
  : Compute Common EOF Patterns from Multiple Sources

## Cross-Validation

Tune hyperparameters

- [`prep_folds()`](https://nick-gauthier.github.io/tidyEOF/reference/prep_folds.md)
  : Prepare contiguous cross-validation folds
- [`prep_cv_folds()`](https://nick-gauthier.github.io/tidyEOF/reference/prep_cv_folds.md)
  : Prepare cross-validation folds for CCA tuning
- [`tune_eof()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_eof.md)
  : Cross-validate EOF truncation for a single field
- [`tune_cca()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_cca.md)
  : Tune CCA hyperparameters via cross-validation
- [`summarize_cv()`](https://nick-gauthier.github.io/tidyEOF/reference/summarize_cv.md)
  : Summarize cross-validation results
- [`summarize_eof_cv()`](https://nick-gauthier.github.io/tidyEOF/reference/summarize_eof_cv.md)
  : Summarize EOF cross-validation results

## Diagnostics

Significance testing and visualization

- [`eigen_test()`](https://nick-gauthier.github.io/tidyEOF/reference/eigen_test.md)
  : Test EOF significance using modified Rule N
- [`get_correlation()`](https://nick-gauthier.github.io/tidyEOF/reference/get_correlation.md)
  : Calculate teleconnections between a patterns object and another
  layer
- [`get_fdr()`](https://nick-gauthier.github.io/tidyEOF/reference/get_fdr.md)
  : Calculate FDR-corrected significance contours for teleconnections

## S3 Methods

Print, plot, predict, and subset methods for tidyeof objects

- [`plot(`*`<patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/plot.patterns.md)
  : Plot method for patterns objects
- [`plot(`*`<coupled_patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/plot.coupled_patterns.md)
  : Plot coupled patterns diagnostics
- [`plot(`*`<common_patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/plot.common_patterns.md)
  : Plot method for common_patterns objects
- [`predict(`*`<coupled_patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/predict.coupled_patterns.md)
  : Predict Method for Coupled Patterns
- [`print(`*`<patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/print.patterns.md)
  : Print method for patterns objects
- [`print(`*`<coupled_patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/print.coupled_patterns.md)
  : Print method for coupled_patterns
- [`print(`*`<cv_folds>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/print.cv_folds.md)
  : Print method for cv_folds objects
- [`screeplot(`*`<patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/screeplot.patterns.md)
  : Scree plot for EOF patterns
- [`summary(`*`<coupled_patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/summary.coupled_patterns.md)
  : Summary method for coupled_patterns
- [`` `[`( ``*`<patterns>`*`)`](https://nick-gauthier.github.io/tidyEOF/reference/sub-.patterns.md)
  : Extract subset of PCs from a patterns object
