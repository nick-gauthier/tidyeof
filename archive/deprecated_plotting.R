# Deprecated plotting functions
# Archived from R/plot_patterns.R
# These are no longer exported but kept for reference

#' Plot EOF spatial patterns (deprecated)
#'
#' @param patterns A patterns object
#' @param scaled Logical, whether to show correlations (TRUE) or raw loadings (FALSE)
#' @param rawdata Optional raw data for correlation calculation when scaled = TRUE
#' @return A ggplot object
#' @seealso \code{\link{plot.patterns}} for the recommended plotting interface
#'
#' @details
#' \strong{This function is deprecated.} Use \code{plot(patterns)} or
#' \code{plot(patterns, type = "eofs")} instead.
plot_eofs <- function(patterns, scaled = FALSE, rawdata = NULL){
  .Deprecated("plot.patterns", package = "tidyeof",
              msg = "plot_eofs() is deprecated. Use plot(patterns) or plot(patterns, type = 'eofs') instead.")

  # Call the new method
  plot.patterns(patterns, type = "eofs", scaled = scaled, rawdata = rawdata)
}

#' Plot PC amplitude time series (deprecated)
#'
#' @param patterns A patterns object
#' @param scale How to scale amplitudes: "standardized" (default), "variance", or "raw"
#' @param scale_y Y-axis scaling: "fixed" (default) or "free"
#' @param events Optional vector of dates to mark with vertical lines
#' @return A ggplot object
#' @seealso \code{\link{plot.patterns}} for the recommended plotting interface
#'
#' @details
#' \strong{This function is deprecated.} Use \code{plot(patterns, type = "amplitudes")}
#' with additional parameters instead.
plot_amps <- function(patterns,
                     scale = c("standardized", "variance", "raw"),
                     scale_y = c("fixed", "free"),
                     events = NULL) {

  .Deprecated("plot.patterns", package = "tidyeof",
              msg = "plot_amps() is deprecated. Use plot(patterns, type = 'amplitudes', scale = '...', scale_y = '...') instead.")

  scale <- match.arg(scale)
  scale_y <- match.arg(scale_y)

  # Call the new method
  plot.patterns(patterns, type = "amplitudes", scale = scale, scale_y = scale_y, events = events)
}

#' Scree plot (deprecated)
#'
#' @param dat Data object
#' @param k Number of components
#' @param kmax Maximum components to show
#' @param scale Scaling option
#' @param monthly Monthly option
#' @param weight Weighting option
#' @details
#' \strong{This function is deprecated.} Use \code{screeplot(pat)} instead
#' where pat is the result of \code{patterns()}.
plot_scree <- function(dat, k = NULL, kmax = 10, scale = FALSE, monthly = FALSE, weight = TRUE) {
  .Deprecated("screeplot.patterns", package = "tidyeof",
              msg = "plot_scree() is deprecated. Use screeplot(patterns(dat, ...)) instead.")

  # For backward compatibility, try to create patterns and call new method
  tryCatch({
    pat <- patterns(dat, k = if(is.null(k)) 5 else k, scale = scale, monthly = monthly, weight = weight)
    screeplot.patterns(pat, k = k, kmax = kmax)
  }, error = function(e) {
    stop("plot_scree() is deprecated and could not convert to new interface. Use screeplot(patterns(dat, ...)) instead.")
  })
}

#' Calculate significance of EOF patterns using North's Rule of Thumb
#'
#' @param patterns A patterns object
#' @param confidence Confidence level (default 0.95)
#' @return A stars object with significance regions
#' @keywords internal
calculate_significance <- function(patterns, confidence = 0.95) {
  # Get key dimensions
  n_times <- nrow(patterns$amplitudes)
  eigenvals <- patterns$eigenvalues$eigenvalues

  # Calculate sampling error for eigenvalues using North's rule
  # Delta lambda = lambda * sqrt(2/n)
  sampling_errors <- eigenvals * sqrt(2/n_times)

  # Test if eigenvalues are distinct
  # If error bars overlap, patterns may be mixed
  is_distinct <- logical(length(eigenvals))
  for(i in seq_along(eigenvals)) {
    if(i < length(eigenvals)) {
      # Check if error bars overlap with next eigenvalue
      delta_i <- sampling_errors[i]
      delta_j <- sampling_errors[i + 1]
      separation <- eigenvals[i] - eigenvals[i + 1]
      is_distinct[i] <- separation > (delta_i + delta_j)
    } else {
      # Last eigenvalue is compared with zero
      is_distinct[i] <- eigenvals[i] > sampling_errors[i]
    }
  }

  # Create significance field for each EOF
  # Use same dimensions as original EOFs
  signif_array <- patterns$eofs[[1]] * 0  # Initialize with zeros

  for(i in seq_len(patterns$k)) {
    if(is_distinct[i]) {
      # If eigenvalue is distinct, mark areas where EOF loading
      # is above threshold as significant
      pattern <- abs(patterns$eofs[[1]][,,i])
      # Use standard error to set threshold
      threshold <- sqrt(1/n_times)  # Simple approximation for significance
      signif_array[,,i] <- pattern > threshold
    }
  }

  # Convert to stars object
  result <- patterns$eofs
  result[[1]] <- signif_array
  names(result) <- "significance"

  # Store significance info in attributes
  attr(result, "is_distinct") <- is_distinct
  attr(result, "sampling_errors") <- sampling_errors
  attr(result, "n_effective") <- n_times

  result
}

#' Print significance test results
#'
#' @param signif_obj Significance object from calculate_significance
print_significance <- function(signif_obj) {
  is_distinct <- attr(signif_obj, "is_distinct")
  sampling_errors <- attr(signif_obj, "sampling_errors")
  n_eff <- attr(signif_obj, "n_effective")

  cat("\nEOF Significance Testing (North's Rule of Thumb)\n")
  cat("Number of effective samples:", n_eff, "\n\n")

  for(i in seq_along(is_distinct)) {
    cat(glue::glue("EOF {i}: {ifelse(is_distinct[i], 'Distinct', 'May be degenerate')} (\\u00B1{round(sampling_errors[i], 3)})\n"))
  }
}
