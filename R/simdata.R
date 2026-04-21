#' Simulate paired generation-structured data
#'
#' `simdata()` generates synthetic paired Gaussian data generation by generation.
#' The output is designed to be passed directly to [treeangle()].
#'
#' The simulated values should be interpreted as paired observations or paired
#' increments, not as a raw nested tree object.
#'
#' Two decay rules are available:
#' - `"power"` uses `rho^i`, matching the manuscript simulation setup;
#' - `"linear"` reproduces the old script logic more closely and uses
#'   `rho * (1 - (i - 1) / (generations - 1))` when `generations > 1`.
#'
#' @param generations Positive integer. Number of generations.
#' @param branching Positive integer. Branching factor. `2` corresponds to a
#' binary tree.
#' @param observations_per_node Positive integer. Number of paired observations
#' contributed by each node in each generation.
#' @param mu Numeric scalar or length-two numeric vector. Mean of the paired
#' Gaussian observations. If a scalar is supplied, it is recycled to both
#' coordinates.
#' @param sigma2 Positive numeric scalar or length-two numeric vector. Marginal
#' variances of the paired Gaussian observations. If a scalar is supplied, it is
#' recycled to both coordinates.
#' @param rho Numeric scalar in `[0, 1)`. Base correlation parameter.
#' @param decay Character string. One of `"power"` or `"linear"`.
#' @param format Character string. One of `"list"`, `"matrix"`, or `"xy"`.
#' @param seed Optional integer seed.
#'
#' @return A `treeangle_simulation` object with components:
#' \itemize{
#'   \item `data`: simulated data in the requested format;
#'   \item `data_by_generation`: list of two-column matrices, one per generation;
#'   \item `correlation_by_generation`: correlation used in each generation;
#'   \item `n_by_generation`: number of observations in each generation;
#'   \item the simulation settings used to generate the data.
#' }
#'
#' @seealso [treeangle()]
#'
#' @examples
#' sim1 <- simdata(seed = 1)
#' print(sim1)
#' treeangle(sim1, alpha = 0.95, method = "mean")
#'
#' sim2 <- simdata(decay = "linear", format = "matrix", seed = 1)
#' head(sim2$data)
#' treeangle(sim2, normalize = TRUE)
#'
#' @export
simdata <- function(generations = 6,
                    branching = 2,
                    observations_per_node = 1,
                    mu = 7,
                    sigma2 = 0.02,
                    rho = 0.5,
                    decay = c("power", "linear"),
                    format = c("list", "matrix", "xy"),
                    seed = NULL) {
  call <- match.call()
  decay <- match.arg(decay)
  format <- match.arg(format)

  generations <- validate_positive_integer_scalar(generations, "generations")
  branching <- validate_positive_integer_scalar(branching, "branching")
  observations_per_node <- validate_positive_integer_scalar(
    observations_per_node,
    "observations_per_node"
  )

  mu <- expand_length_two(mu, "mu")
  sigma2 <- expand_length_two(sigma2, "sigma2")

  if (any(!is.finite(mu))) {
    stop("'mu' must contain finite numeric values.", call. = FALSE)
  }

  if (any(!is.finite(sigma2)) || any(sigma2 <= 0)) {
    stop("'sigma2' must contain positive finite numeric values.", call. = FALSE)
  }

  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || rho < 0 || rho >= 1) {
    stop("'rho' must be a single numeric value in [0, 1).", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != floor(seed)) {
      stop("'seed' must be NULL or a single integer value.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  generation_index <- seq_len(generations)
  n_by_generation <- observations_per_node * branching^generation_index

  if (any(n_by_generation > .Machine$integer.max)) {
    stop(
      "The requested simulation size is too large for integer indexing. ",
      "Please reduce 'generations', 'branching', or 'observations_per_node'.",
      call. = FALSE
    )
  }

  n_by_generation <- as.integer(n_by_generation)
  correlation_by_generation <- correlation_sequence(
    generations = generations,
    rho = rho,
    decay = decay
  )

  data_by_generation <- lapply(seq_len(generations), function(i) {
    sigma_mat <- covariance_matrix(
      sigma2 = sigma2,
      correlation = correlation_by_generation[i]
    )

    sim_i <- rmvnorm_chol(
      n = n_by_generation[i],
      mean = mu,
      sigma = sigma_mat
    )

    colnames(sim_i) <- c("x", "y")
    sim_i
  })

  data_out <- switch(
    format,
    list = data_by_generation,
    matrix = bind_generation_matrix(data_by_generation),
    xy = do.call(rbind, data_by_generation)
  )

  if (format == "xy") {
    colnames(data_out) <- c("x", "y")
  }

  out <- list(
    data = data_out,
    data_by_generation = data_by_generation,
    format = format,
    generations = generations,
    branching = branching,
    observations_per_node = observations_per_node,
    mu = mu,
    sigma2 = sigma2,
    rho = rho,
    decay = decay,
    correlation_by_generation = correlation_by_generation,
    n_by_generation = n_by_generation,
    total_observations = sum(n_by_generation),
    call = call
  )

  class(out) <- "treeangle_simulation"
  out
}

#' @export
#' @noRd
print.treeangle_simulation <- function(x, ...) {
  cat("<treeangle_simulation>\n")
  cat("  format                : ", x$format, "\n", sep = "")
  cat("  generations           : ", x$generations, "\n", sep = "")
  cat("  branching             : ", x$branching, "\n", sep = "")
  cat("  observations_per_node : ", x$observations_per_node, "\n", sep = "")
  cat("  rho                   : ", formatC(x$rho, digits = 4, format = "f"), "\n", sep = "")
  cat("  decay                 : ", x$decay, "\n", sep = "")
  cat("  total_observations    : ", x$total_observations, "\n", sep = "")
  invisible(x)
}

validate_positive_integer_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x != floor(x)) {
    stop(
      "'",
      arg_name,
      "' must be a single positive integer.",
      call. = FALSE
    )
  }

  as.integer(x)
}

expand_length_two <- function(x, arg_name) {
  if (!is.numeric(x) || anyNA(x)) {
    stop(
      "'",
      arg_name,
      "' must be numeric and must not contain missing values.",
      call. = FALSE
    )
  }

  if (length(x) == 1L) {
    return(rep(as.numeric(x), 2L))
  }

  if (length(x) == 2L) {
    return(as.numeric(x))
  }

  stop(
    "'",
    arg_name,
    "' must have length 1 or 2.",
    call. = FALSE
  )
}

correlation_sequence <- function(generations, rho, decay) {
  index <- seq_len(generations)

  if (decay == "power") {
    return(rho^index)
  }

  if (generations == 1L) {
    return(rho)
  }

  rho * (1 - (index - 1) / (generations - 1))
}

covariance_matrix <- function(sigma2, correlation) {
  cov12 <- correlation * sqrt(sigma2[1] * sigma2[2])

  matrix(
    c(
      sigma2[1], cov12,
      cov12, sigma2[2]
    ),
    nrow = 2,
    byrow = TRUE
  )
}

bind_generation_matrix <- function(data_list) {
  out <- do.call(
    rbind,
    lapply(seq_along(data_list), function(i) {
      data.frame(
        x = data_list[[i]][, 1],
        y = data_list[[i]][, 2],
        generation = i
      )
    })
  )

  rownames(out) <- NULL
  out
}

rmvnorm_chol <- function(n, mean, sigma) {
  sigma <- (sigma + t(sigma)) / 2
  z <- matrix(stats::rnorm(n * length(mean)), nrow = n, ncol = length(mean))
  chol_sigma <- chol(sigma)
  out <- z %*% chol_sigma
  out <- sweep(out, 2, mean, FUN = "+")
  storage.mode(out) <- "double"
  out
}