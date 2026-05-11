#' Simulate paired generation-structured data
#'
#' `simdata()` generates synthetic paired Gaussian data generation by generation.
#' It returns a plain three-column data frame with columns:
#' `generation`, `x`, and `y`.
#'
#' The output is designed to be passed directly to [treeangle()].
#'
#' The simulated values should be interpreted as paired observations or paired
#' increments, not as a raw nested tree object.
#'
#' Two decay rules are available:
#'
#' - `"power"` uses `rho^i`, matching the manuscript simulation setup;
#' - `"linear"` uses the manuscript-style linear rule
#'   `rho * (1 - (i - 1) / generations)`.
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
#' @param format Deprecated. Kept only for backward compatibility. Ignored.
#' `simdata()` always returns a data frame.
#' @param seed Optional integer seed.
#'
#' @return A data frame with three columns:
#' \itemize{
#'   \item `generation`: generation index;
#'   \item `x`: first coordinate;
#'   \item `y`: second coordinate.
#' }
#'
#' The returned data frame also contains optional attributes
#' `n_by_generation`, `correlation_by_generation`, and `simulation_settings`.
#'
#' @seealso [treeangle()]
#'
#' @examples
#' dat <- simdata(seed = 1)
#' head(dat)
#' treeangle(dat, alpha = 0.95, method = "mean")
#'
#' @export
simdata <- function(generations = 6,
                    branching = 2,
                    observations_per_node = 1,
                    mu = 7,
                    sigma2 = 0.02,
                    rho = 0.5,
                    decay = c("power", "linear"),
                    format = c("data.frame", "list", "matrix", "xy"),
                    seed = NULL) {
  call <- match.call()
  decay <- match.arg(decay)

  # Kept only for backward compatibility. It no longer controls the output type.
  requested_format <- match.arg(format)

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

  if (
    !is.numeric(rho) ||
      length(rho) != 1L ||
      is.na(rho) ||
      !is.finite(rho) ||
      rho < 0 ||
      rho >= 1
  ) {
    stop("'rho' must be a single numeric value in [0, 1).", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (
      !is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed != floor(seed)
    ) {
      stop("'seed' must be NULL or a single integer value.", call. = FALSE)
    }

    set.seed(as.integer(seed))
  }

  generation_index <- seq_len(generations)
  n_by_generation <- observations_per_node * branching^generation_index

  if (
    any(!is.finite(n_by_generation)) ||
      any(n_by_generation > .Machine$integer.max)
  ) {
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

  names(data_by_generation) <- paste0("generation_", seq_len(generations))

  out <- do.call(
    rbind,
    lapply(seq_along(data_by_generation), function(i) {
      z <- data_by_generation[[i]]

      data.frame(
        generation = rep.int(i, nrow(z)),
        x = z[, 1],
        y = z[, 2]
      )
    })
  )

  row.names(out) <- NULL
  out$generation <- as.integer(out$generation)

  attr(out, "n_by_generation") <- n_by_generation
  attr(out, "correlation_by_generation") <- correlation_by_generation
  attr(out, "simulation_settings") <- list(
    generations = generations,
    branching = branching,
    observations_per_node = observations_per_node,
    mu = mu,
    sigma2 = sigma2,
    rho = rho,
    decay = decay,
    requested_format = requested_format,
    total_observations = sum(n_by_generation),
    call = call
  )

  out
}

validate_positive_integer_scalar <- function(x, arg_name) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0 ||
      x != floor(x)
  ) {
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

  rho * (1 - (index - 1) / generations)
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

rmvnorm_chol <- function(n, mean, sigma) {
  sigma <- (sigma + t(sigma)) / 2

  z <- matrix(
    stats::rnorm(n * length(mean)),
    nrow = n,
    ncol = length(mean)
  )

  chol_sigma <- chol(sigma)
  out <- z %*% chol_sigma
  out <- sweep(out, 2, mean, FUN = "+")

  storage.mode(out) <- "double"
  out
}