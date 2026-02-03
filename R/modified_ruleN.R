#' Test EOF significance using modified Rule N
#'
#' Tests whether the k-th eigenvalue is significantly different from noise
#' using a modified Rule N approach based on Tracy-Widom distribution.
#'
#' @param lambdas Vector of eigenvalues from PCA
#' @param k Index of eigenvalue to test
#' @param M Number of spatial points (grid cells)
#' @param n Number of time steps
#'
#' @return Logical, TRUE if eigenvalue is significant at p < 0.05
#' @export
eigen_test <- function(lambdas, k, M, n){
  nrank<- min(n - 1, M)

  kstar <- M - k + 1
  nk <- n - k + 1
  mu <- (sqrt(nk - 0.5) + sqrt(kstar - 0.5))^2
  sigma <- sqrt(mu) * (1 / sqrt(nk - 0.5) + 1 / sqrt(kstar - 0.5)) ^ (1/3)
  alpha <- 46.4
  beta <- (0.186 * sigma) / max(nk, kstar)
  zeta <- (mu - 9.85 * sigma) / max(nk, kstar)
  lambda_star <- lambdas[k] / ((1 / (nrank - k + 1)) * sum(lambdas[k:nrank]))


  (1- pgamma(((lambda_star - zeta) / beta), 46.4)) < 0.05
}

#test2 <- c(3.532, 1.985, .344, .074, .038, .027)


#map_lgl(1:6, ~eigen_test(test2, k = ., M = 6, n = 31))


#map_lgl(1:36, ~eigen_test(prism_eigs$eigenvalues, k = ., M = 17783, n = 36))

#map_lgl(1:275, ~eigen_test(sim_eigs$eigenvalues, k = ., M = 276, n = 1001))
