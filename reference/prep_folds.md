# Prepare contiguous cross-validation folds

Divides time steps into k contiguous folds for temporal
cross-validation. Folds are kept contiguous to preserve temporal
structure.

## Usage

``` r
prep_folds(times, kfolds = 5)
```

## Arguments

- times:

  Vector of time values to split

- kfolds:

  Number of folds (default 5)

## Value

List of k vectors, each containing the time values for that fold

## See also

\[prep_cv_folds()\] for preparing complete CV folds with pre-computed
patterns

\[tune_cca()\] for hyperparameter grid search
