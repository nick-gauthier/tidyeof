# Archived Legacy Code

These files contain legacy code from earlier versions of tidyEOF that used different
approaches (GAM-based coupling, remote package for EOT, etc.). Kept for reference
when replicating old analyses.

## Files

- `predict_patterns.R` - Old GAM-based prediction using modelr::add_predictions
- `validate.R` - Legacy cross-validation error metrics (hardcoded for specific datasets)

## See Also

- `R/legacy_wrappers.R` - Commented-out legacy wrapper functions
- `R/prep_cv.R` - Contains commented-out `prep_eot()` that used remote package
