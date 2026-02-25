# Changelog

## tidyeof 0.1.0

Initial release.

### Core Features

- [`patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/patterns.md)
  for EOF extraction from `stars` spatiotemporal data, with area
  weighting, standardization, varimax rotation, and monthly
  anomalization
- [`couple()`](https://nick-gauthier.github.io/tidyEOF/reference/couple.md)
  and [`predict()`](https://rdrr.io/r/stats/predict.html) for CCA-based
  statistical downscaling
- [`common_patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/common_patterns.md)
  for joint EOF decomposition across multiple datasets
- [`tune_eof()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_eof.md)
  and
  [`tune_cca()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_cca.md)
  for cross-validated hyperparameter selection with contiguous temporal
  folds
- Scree plots with North et al. (1982) error bars and modified Rule N
  significance testing
- Teleconnection maps with FDR-corrected significance contours
- Full integration with `stars`, `sf`, and the tidyverse

### Architecture

- S3 classes: `patterns`, `coupled_patterns`, `common_patterns`,
  `cv_folds`
- Plot, print, summary, and predict methods for all major classes
- Automatic IRLBA truncated SVD for large datasets
- Unit preservation throughout the pipeline via the `units` package
