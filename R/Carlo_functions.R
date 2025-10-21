#' @useDynLib StratMonteCarlo, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils tail
NULL

# =============================================================
# MONTE CARLO INTEGRATION FUNCTION
# =============================================================

#' Monte Carlo Integration with Optional Importance Sampling
#'
#' @description
#' Performs Monte Carlo integration over 1D or multi-dimensional domains.
#' Supports function inputs or expression strings (evaluated via tinyexpr),
#' optional adaptive partitioning, and importance sampling using different distributions.
#'
#' @param f Function or character expression representing the integrand.
#' @param lower Numeric vector of lower bounds.
#' @param upper Numeric vector of upper bounds.
#' @param n_samples Integer; number of Monte Carlo samples.
#' @param partition Logical; whether to use adaptive stratified partitioning.
#' @param dim Integer; integration dimension. If NULL, inferred from length(lower).
#' @param step Integer; number of samples per recorded estimate (for plotting convergence).
#' @param importance_sampling Logical; whether to use importance sampling.
#' @param is_distribution Character; importance sampling distribution, one of "normal", "beta", "exponential", "mixture_normal".
#' @param is_params Optional list; user-supplied importance sampling parameters.
#' @param expr_vars Optional character vector; variable names for expression input.
#'
#' @return An object of class \code{"mc_result"} containing:
#'   \item{estimate}{Final Monte Carlo estimate}
#'   \item{estimates}{Vector of estimates at each step for convergence plotting}
#'   \item{importance_sampling}{Logical; whether importance sampling was used}
#'   \item{is_distribution}{Distribution used for importance sampling (if any)}
#'   \item{params_used}{Optimized or user-supplied IS parameters (if importance_sampling = TRUE)}
#'   \item{expr}{Expression used (if f is character)}
#'   \item{n_partitions}{Number of partitions used (if partitioning enabled)}
#'
#' @examples
#' f <- function(x) x^2
#' montecarlo_integrate(f, lower = 0, upper = 1, n_samples = 1e4)
#'
#' f_expr <- "x1^2 + x2^2"
#' montecarlo_integrate(f_expr, lower = c(0,0), upper = c(1,1), n_samples = 1e4)
#'
#' @export
montecarlo_integrate <- function(f, lower, upper, n_samples,
                                 partition = FALSE, dim = NULL,
                                 step = 100, importance_sampling = FALSE,
                                 is_distribution = c("normal", "beta", "exponential", "mixture_normal"),
                                 is_params = NULL, expr_vars = NULL) {

  is_distribution <- match.arg(is_distribution)
  if (is.null(dim)) dim <- length(lower)

  if (!is.function(f) && !is.character(f))
    stop("'f' must be a function or character expression.")

  if (is.character(f) && is.null(expr_vars))
    expr_vars <- paste0("x", seq_len(dim))

  # Convert character expression to function for partitioning
  f_for_partition <- NULL
  if (partition && is.character(f)) {
    if (dim == 1) {
      expr_clean <- gsub("x1", "x", f)
      f_for_partition <- eval(parse(text = paste0("function(x) { ", expr_clean, " }")))
    } else {
      var_assigns <- paste0(expr_vars, " <- x[", seq_len(dim), "]", collapse = "; ")
      func_text <- paste0("function(x) { ", var_assigns, "; ", f, " }")
      f_for_partition <- eval(parse(text = func_text))
    }
  } else if (partition && is.function(f)) {
    f_for_partition <- f
  }

  result <- montecarlo_integrate_cpp(
    f_input = if (!is.null(f_for_partition)) f_for_partition else if (is.function(f)) f else NULL,
    expr = if (is.character(f) && !partition) f else "",
    vars = expr_vars,
    lower = lower,
    upper = upper,
    n_samples = n_samples,
    partition = partition,
    dim = dim,
    importance_sampling = importance_sampling,
    is_distribution = is_distribution,
    is_params = is_params,
    step = step
  )

  class(result) <- "mc_result"
  result
}

# =============================================================
# PRINT METHOD
# =============================================================

#' Print Monte Carlo Result
#'
#' @param x Object of class \code{"mc_result"}
#' @param ... Ignored
#' @export
print.mc_result <- function(x, ...) {
  cat("Monte Carlo Integration Result\n")
  cat("================================\n")
  cat("Final estimate:", formatC(x$estimate, digits = 6, format = "g"), "\n")
  cat("Number of recorded estimates:", length(x$estimates), "\n")
  if (!is.null(x$n_partitions) && x$n_partitions > 1) cat("Partitions:", x$n_partitions, "\n")
  if (x$importance_sampling) {
    cat("Importance Sampling:", x$is_distribution, "\n")
    if (!is.null(x$params_used)) {
      cat("Parameters:\n")
      print(x$params_used)
    }
  }
  invisible(x)
}

# =============================================================
# PLOT FUNCTION
# =============================================================

#' Plot Monte Carlo Estimates Convergence
#'
#' @param ... Named mc_result objects
#' @param step Number of samples per estimate (for x-axis scaling)
#' @param true_value Optional numeric; true integral value for reference line
#' @param title_size Title text size
#' @param axis_title_size Axis title text size
#' @param axis_text_size Axis label text size
#' @param legend_text_size Legend text size
#' @return ggplot2 object
#' @export
mc_plot <- function(..., step = 100, true_value = NULL,
                    title_size = 14, axis_title_size = 12,
                    axis_text_size = 10, legend_text_size = 10) {

  # Avoid NSE notes from R CMD check
  Samples <- Estimate <- Method <- NULL

  args <- list(...)
  arg_names <- names(args)
  if (is.null(arg_names)) arg_names <- paste0("Method", seq_along(args))
  arg_names[arg_names == ""] <- paste0("Method", which(arg_names == ""))

  df_list <- lapply(seq_along(args), function(i) {
    obj <- args[[i]]
    est_vec <- if (inherits(obj, "mc_result")) obj$estimates else if (is.numeric(obj)) obj else stop("Invalid object")
    data.frame(
      Samples = seq_len(length(est_vec)) * step,
      Estimate = as.numeric(est_vec),
      Method = arg_names[i],
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Samples, y = Estimate, color = Method)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Monte Carlo Estimates Convergence",
                  x = "Number of Samples", y = "Estimate") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold"),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_brewer(palette = "Set1")

  if (!is.null(true_value)) {
    p <- p +
      ggplot2::geom_hline(yintercept = true_value, linetype = "dashed", color = "black", linewidth = 0.8) +
      ggplot2::annotate("text",
                        x = max(df$Samples) * 0.95,
                        y = true_value,
                        label = paste("True =", round(true_value, 3)),
                        vjust = -0.5, hjust = 1, color = "black", size = 4)
  }

  p
}

# =============================================================
# PRINT FINAL ESTIMATE ONLY
# =============================================================

#' Print Only Final Monte Carlo Estimate
#'
#' @param x mc_result or numeric vector
#' @export
mc_print <- function(x) {
  val <- if (inherits(x, "mc_result")) x$estimate else if (is.numeric(x)) tail(x, 1) else stop("Unsupported object")
  cat("Monte Carlo estimate (final):", formatC(val, digits = 6, format = "g"), "\n")
  invisible(val)
}

# =============================================================
# PREVIEW IMPORTANCE SAMPLING PARAMETERS
# =============================================================

#' Preview Optimized Importance Sampling Parameters
#'
#' @param f Function or expression
#' @param lower Numeric vector; lower bounds
#' @param upper Numeric vector; upper bounds
#' @param dim Integer; dimension (default inferred)
#' @param distribution Character; one of "normal", "beta", "exponential", "mixture_normal"
#' @param n_pilot Integer; number of pilot samples
#' @export
preview_is_params <- function(f, lower, upper, dim = NULL,
                              distribution = c("normal", "beta", "exponential", "mixture_normal"),
                              n_pilot = 1000) {
  distribution <- match.arg(distribution)
  if (is.null(dim)) dim <- length(lower)
  if (length(lower) == 1) lower <- rep(lower, dim)
  if (length(upper) == 1) upper <- rep(upper, dim)

  cat("Parameters are optimized automatically during integration.\n")
  cat("To inspect them, check the $params_used field of the mc_result object.\n\n")
  cat("Example:\n")
  cat("  res <- montecarlo_integrate(f, lower, upper, n_samples = 1e4, importance_sampling = TRUE)\n")
  cat("  res$params_used\n")
  invisible(NULL)
}


#' Extract Importance Sampling Parameters
#'
#' @param mc_res An object of class "mc_result"
#' @return A list of importance sampling parameters used for each dimension
#' @export
get_is_params <- function(mc_res) {
  if (!inherits(mc_res, "mc_result")) stop("Input must be an mc_result object.")
  if (!mc_res$importance_sampling) stop("This result was not generated using importance sampling.")

  return(mc_res$params_used)
}
