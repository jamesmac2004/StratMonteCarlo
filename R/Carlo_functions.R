#' Wrapper for Monte Carlo integration (vector of cumulative estimates)
#'
#' @param f Function to integrate
#' @param lower Lower bound(s)
#' @param upper Upper bound(s)
#' @param n_samples Number of samples
#' @param Partition Use stratified partitioning (logical, not yet used in this version)
#' @param dim Dimension of integration
#' @param step Interval at which to record cumulative estimates
#' @return An object of class "mc_result" containing:
#'   - estimates: numeric vector of cumulative estimates
#'   - final: last Monte Carlo estimate
#' @export
montecarlo_integrate <- function(f, lower, upper, n_samples, Partition = FALSE, dim = 1, step = 100) {
  stopifnot(is.function(f))
  stopifnot(is.numeric(lower))
  stopifnot(is.numeric(upper))
  stopifnot(is.numeric(n_samples) && n_samples > 0)
  stopifnot(is.logical(Partition))

  if (length(lower) == 1) lower <- rep(lower, dim)
  if (length(upper) == 1) upper <- rep(upper, dim)

  volume <- prod(upper - lower)
  n_steps <- n_samples %/% step
  estimates <- numeric(n_steps)

  running_sum <- 0
  x <- numeric(dim)

  for (i in seq_len(n_samples)) {
    for (d in seq_len(dim)) {
      x[d] <- lower[d] + (upper[d] - lower[d]) * runif(1)
    }

    val <- if (dim == 1) f(x[1]) else do.call(f, as.list(x))
    running_sum <- running_sum + val

    if (i %% step == 0) {
      estimates[i %/% step] <- volume * running_sum / i
    }
  }

  # 👇 return as structured object, not a raw vector
  result <- list(
    estimates = estimates,
    final = tail(estimates, 1)
  )
  class(result) <- "mc_result"
  return(result)
}

#' Print method for Monte Carlo results
#' @param x An object of class "mc_result"
#' @param ... Ignored
#' @export
print.mc_result <- function(x, ...) {
  cat("Monte Carlo estimate (final):", x$final, "\n")
}


#' Plot precomputed Monte Carlo estimates
#'
#' @param ... Named mc_result objects
#' @param step Number of samples per estimate (for x-axis scaling)
#' @param true_value Optional, true integral value
#' @return ggplot of MC convergence
#' @export
mc_plot <- function(..., step = 1, true_value = NULL) {
  library(ggplot2)
  estimates_list <- list(...)

  # Extract estimates if objects are mc_result
  estimates_list <- lapply(estimates_list, function(obj) {
    if (inherits(obj, "mc_result")) obj$estimates else obj
  })

  lengths_vec <- sapply(estimates_list, length)
  if(length(unique(lengths_vec)) != 1) stop("All inputs must have the same length")

  df <- do.call(rbind, lapply(names(estimates_list), function(nm) {
    data.frame(
      Samples = seq(step, step * lengths_vec[1], by = step),
      Estimate = estimates_list[[nm]],
      Method = nm
    )
  }))

  p <- ggplot(df, aes(x = Samples, y = Estimate, color = Method)) +
    geom_line(size = 0.5) +
    geom_point(size = 0.8, alpha = 0.3) +
    theme_minimal() +
    labs(title = "Monte Carlo Estimates Convergence",
         x = "Number of Samples",
         y = "Estimate") +
    theme(text = element_text(size = 14),
          legend.position = "bottom",
          legend.title = element_blank()) +
    scale_color_brewer(palette = "Set1")

  if(!is.null(true_value)) {
    p <- p + geom_hline(yintercept = true_value, linetype = "dashed", color = "black", size = 0.6) +
      annotate("text", x = max(df$Samples)*0.95, y = true_value,
               label = paste("True value =", signif(true_value, 5)),
               vjust = -0.5, hjust = 1, color = "black", size = 4)
  }

  return(p)
}


#' Print only the final Monte Carlo estimate
#'
#' Accepts an object of class "mc_result", a list containing an `estimates`
#' numeric vector, or a raw numeric vector. Prints only the final estimate.
#'
#' @param x mc_result (or list with $estimates) or numeric vector
#' @export
mc_print <- function(x) {
  # 1) mc_result with $final
  if (is.list(x) && inherits(x, "mc_result") && !is.null(x$final)) {
    val <- x$final
  } else if (is.list(x) && !is.null(x$estimates) && is.numeric(x$estimates)) {
    # 2) generic list that contains estimates
    val <- tail(x$estimates, 1)
  } else if (is.numeric(x)) {
    # 3) numeric vector (or scalar)
    val <- tail(x, 1)
  } else {
    stop("mc_print: unsupported object. Provide an mc_result, a list with $estimates, or a numeric vector.")
  }

  cat("Monte Carlo estimate (final):", formatC(val, digits = 6, format = "g"), "\n")
  invisible(val)
}
