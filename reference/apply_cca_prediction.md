# Apply CCA Prediction Transform

Internal function that applies the CCA transformation to predict
response amplitudes from predictor amplitudes.

## Usage

``` r
apply_cca_prediction(new_amplitudes, cca_result, k)
```

## Arguments

- new_amplitudes:

  Tibble with time and predictor PC amplitudes

- cca_result:

  CCA result object from cancor()

- k:

  Number of CCA modes to use

## Value

Tibble with time and predicted response PC amplitudes
