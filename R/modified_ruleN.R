#' Test EOF significance using modified Rule N
#'
#' Tests whether the k-th eigenvalue is significantly different from noise
#' using a modified Rule N approach based on Tracy-Widom distribution.
#'
#' Based on the gamma approximation to the Tracy-Widom distribution described
#' in Cheng & Wallace (1993) and implemented following Overland & Preisendorfer
#' (1982). Constants (shape = 46.4, scale factor = 0.186, location = 9.85)
#' derive from fitting the gamma CDF to the Tracy-Widom Type 1 distribution.
#'
#' @param lambdas Vector of eigenvalues from PCA
#' @param k Index of eigenvalue to test
#' @param M Number of spatial points (grid cells)
#' @param n Number of time steps
#' @param p Significance level (default 0.05)
#'
#' @return Logical, TRUE if eigenvalue is significant at level p
#' @export
eigen_test <- function(lambdas, k, M, n, p = 0.05){
  nrank <- min(n - 1, M)

  kstar <- M - k + 1
  nk <- n - k + 1
  mu <- (sqrt(nk - 0.5) + sqrt(kstar - 0.5))^2
  sigma <- sqrt(mu) * (1 / sqrt(nk - 0.5) + 1 / sqrt(kstar - 0.5)) ^ (1/3)
  shape <- 46.4
  beta <- (0.186 * sigma) / max(nk, kstar)
  zeta <- (mu - 9.85 * sigma) / max(nk, kstar)
  lambda_star <- lambdas[k] / ((1 / (nrank - k + 1)) * sum(lambdas[k:nrank]))

  (1 - pgamma(((lambda_star - zeta) / beta), shape)) < p
}

#' Find Rule N significance cutoff for a patterns object
#'
#' Tests each eigenvalue in sequence and returns the index of the last
#' significant mode. This is the maximum k supported by the data according
#' to the modified Rule N test.
#'
#' @param x A patterns object
#' @param p Significance level (default 0.05)
#'
#' @return Integer: index of last significant eigenvalue, or 0 if none are significant
#' @keywords internal
rule_n_cutoff <- function(x, p = 0.05) {
  lambdas <- x$eigenvalues$eigenvalues
  n_times <- nrow(x$amplitudes)
  n_valid <- length(x$valid_pixels)
  max_test <- min(length(lambdas), n_times - 1)

  significant <- vapply(seq_len(max_test), function(k) {
    eigen_test(lambdas, k = k, M = n_valid, n = n_times, p = p)
  }, logical(1))

  # Last TRUE in contiguous sequence from the start
  if (!significant[1]) return(0L)
  # Find where significance first drops off
  first_nonsig <- which(!significant)[1]
  if (is.na(first_nonsig)) max_test else first_nonsig - 1L
}
