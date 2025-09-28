#' Plot EOF spatial patterns (deprecated)
#'
#' @param patterns A patterns object
#' @param scaled Logical, whether to show correlations (TRUE) or raw loadings (FALSE)
#' @param rawdata Optional raw data for correlation calculation when scaled = TRUE
#' @return A ggplot object
#' @export
#' @seealso \code{\link{plot.patterns}} for the recommended plotting interface
#'
#' @details
#' \strong{This function is deprecated.} Use \code{plot(patterns)} or
#' \code{plot(patterns, type = "eofs")} instead.
plot_eofs <- function(patterns, scaled = FALSE, rawdata = NULL){
  .Deprecated("plot.patterns", package = "tidyEOF",
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
#' @export
#' @seealso \code{\link{plot.patterns}} for the recommended plotting interface
#'
#' @details
#' \strong{This function is deprecated.} Use \code{plot(patterns, type = "amplitudes")}
#' with additional parameters instead.
plot_amps <- function(patterns,
                     scale = c("standardized", "variance", "raw"),
                     scale_y = c("fixed", "free"),
                     events = NULL) {

  .Deprecated("plot.patterns", package = "tidyEOF",
              msg = "plot_amps() is deprecated. Use plot(patterns, type = 'amplitudes', scale = '...', scale_y = '...') instead.")

  scale <- match.arg(scale)
  scale_y <- match.arg(scale_y)

  # Call the new method
  plot.patterns(patterns, type = "amplitudes", scale = scale, scale_y = scale_y, events = events)
}

#' Scree plot for EOF patterns
#'
#' @param x A patterns object from get_patterns()
#' @param k Optional number of components to highlight with vertical line
#' @param kmax Maximum number of components to show (default 10)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' patterns <- get_patterns(data, k = 5)
#' screeplot(patterns)
#' screeplot(patterns, k = 3, kmax = 8)
#' }
screeplot.patterns <- function(x, k = NULL, kmax = 10, ...) {
  x$eigenvalues %>%
    dplyr::mutate(separated = if_else(is.na(lag(low)), TRUE, hi < lag(low)),
           multiplet = as.factor(cumsum(separated))) %>%
    filter(PC <= kmax) %>%
    ggplot2::ggplot(aes(x = PC, y = percent * 100)) +
    ggplot2::geom_linerange(aes(x = PC, ymin = low, ymax = hi)) +
    ggplot2::geom_point(size = 2, aes(color = multiplet)) +
    ggplot2::geom_text(aes(x = PC, y = cumvar_line, label = glue::glue("{round(cumulative * 100, 0)}%")), size = 2.5, vjust = 0) +
    ggplot2::labs(x = "Principal Component", y = "Normalized Eigenvalue") +
    ggplot2::geom_vline(xintercept = k + .5, linetype = 2, color = 'red', alpha = .7) +
    ggplot2::theme_bw() +
    ggplot2::guides(color = 'none') +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::scale_color_brewer(palette = 'Spectral')
}

#' Scree plot (deprecated)
#'
#' @param dat Data object
#' @param k Number of components
#' @param kmax Maximum components to show
#' @param scale Scaling option
#' @param monthly Monthly option
#' @param weight Weighting option
#' @export
#' @details
#' \strong{This function is deprecated.} Use \code{screeplot(patterns)} instead
#' where patterns is the result of \code{get_patterns()}.
plot_scree <- function(dat, k = NULL, kmax = 10, scale = FALSE, monthly = FALSE, weight = TRUE) {
  .Deprecated("screeplot.patterns", package = "tidyEOF",
              msg = "plot_scree() is deprecated. Use screeplot(get_patterns(dat, ...)) instead.")

  # For backward compatibility, try to create patterns and call new method
  tryCatch({
    patterns <- get_patterns(dat, k = if(is.null(k)) 5 else k, scale = scale, monthly = monthly, weight = weight)
    screeplot.patterns(patterns, k = k, kmax = kmax)
  }, error = function(e) {
    stop("plot_scree() is deprecated and could not convert to new interface. Use screeplot(get_patterns(dat, ...)) instead.")
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
#' @export
print_significance <- function(signif_obj) {
  is_distinct <- attr(signif_obj, "is_distinct")
  sampling_errors <- attr(signif_obj, "sampling_errors")
  n_eff <- attr(signif_obj, "n_effective")

  cat("\nEOF Significance Testing (North's Rule of Thumb)\n")
  cat("Number of effective samples:", n_eff, "\n\n")

  for(i in seq_along(is_distinct)) {
    cat(glue("EOF {i}: {ifelse(is_distinct[i], 'Distinct', 'May be degenerate')} (±{sampling_errors[i]:.3f})\\n"))
  }
}
