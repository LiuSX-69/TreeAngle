#' Compute the TreeAngle statistic
#'
#' `treeangle()` computes the geometric TreeAngle statistic from paired
#' two-dimensional data.
#'
#' Each row of the supplied data should represent one paired observation
#' `(x, y)`. For generation-structured input, rows are grouped by generation and
#' the function can optionally apply the normalization step described in the
#' manuscript before computing the angle.
#'
#' Supported input formats are:
#' - a list of generation-specific matrices/data frames, each with two numeric
#'   columns;
#' - a matrix/data frame with columns `x`, `y`, and `generation`;
#' - a direct two-column matrix/data frame.
#'
#' Objects returned by [simdata()] can be passed directly to `treeangle()`.
#'
#' @param data Input data. Supported formats are:
#' \itemize{
#'   \item a list of generation-specific matrices/data frames, each with at least
#'   two numeric columns;
#'   \item a matrix/data frame with at least three columns, where columns 1 and 2
#'   are `x`, `y`, and column 3 is a positive integer generation index;
#'   \item a matrix/data frame with at least two numeric columns, in which case
#'   only columns 1 and 2 are used.
#' }
#' @param alpha Numeric scalar in `(0, 1)`. This is the target coverage
#' proportion. For example, `alpha = 0.95` means that about 95\% of points are
#' enclosed by the estimated angle.
#' @param method Character string specifying the summary rule used to choose the
#' final result. Must be one of `"mean"`, `"min"`, `"median"`, `"equal"`, or
#' `"neigh"`.
#' @param normalize Logical scalar or `NULL`. If `NULL` (default), normalization
#' is applied automatically for generation-structured input and skipped for
#' direct two-column input. Set `TRUE` or `FALSE` to override this behavior.
#' @param target_sd Positive numeric scalar. Target marginal standard deviation
#' used in the normalization step. Used only when normalization is applied.
#' Default is `0.01`.
#' @param tau Non-negative numeric scalar. Generation shift parameter used in the
#' normalization step. Used only when normalization is applied. Default is `0.1`.
#'
#' @return A `treeangle_result` object. Important components are:
#' \itemize{
#'   \item `angle_degrees`: estimated angle in degrees;
#'   \item `tan_angle`: tangent of the estimated angle;
#'   \item `upper_point`, `lower_point`: boundary points;
#'   \item `upper_line`, `lower_line`: slope and intercept of the two boundary lines;
#'   \item `method_table`: a table of the five available summary methods;
#'   \item `result_matrix`: legacy 5 x 2 matrix representation.
#' }
#'
#' Use `print()` for a compact summary and `as.matrix()` to recover the legacy
#' matrix form.
#'
#' @seealso [simdata()]
#'
#' @examples
#' sim_list <- simdata(
#'   generations = 6,
#'   rho = 0.7,
#'   decay = "linear",
#'   format = "list",
#'   seed = 123
#' )
#'
#' res1 <- treeangle(
#'   sim_list,
#'   alpha = 0.95,
#'   method = "mean"
#' )
#' print(res1)
#'
#' sim_matrix <- simdata(
#'   generations = 6,
#'   rho = 0.7,
#'   decay = "power",
#'   format = "matrix",
#'   seed = 123
#' )
#'
#' res2 <- treeangle(
#'   sim_matrix,
#'   alpha = 0.95,
#'   method = "median",
#'   normalize = TRUE
#' )
#' print(res2)
#' head(res2$method_table)
#'
#' @export
treeangle <- function(data,
                      alpha = 0.95,
                      method = "mean",
                      normalize = NULL,
                      target_sd = 0.01,
                      tau = 0.1) {
  call <- match.call()
  valid_methods <- c("mean", "min", "median", "equal", "neigh")

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop(
      "'alpha' must be a single numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.character(method) || length(method) != 1L || !method %in% valid_methods) {
    stop(
      "'method' must be one of: ",
      paste(valid_methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(normalize) &&
      !(is.logical(normalize) && length(normalize) == 1L && !is.na(normalize))) {
    stop("'normalize' must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  if (inherits(data, "treeangle_simulation")) {
    data <- data$data
  }

  input_type <- detect_input_type(data)
  normalize_flag <- resolve_normalize_flag(
    normalize = normalize,
    input_type = input_type
  )

  if (normalize_flag) {
    if (!is.numeric(target_sd) || length(target_sd) != 1L || is.na(target_sd) ||
        target_sd <= 0) {
      stop(
        "'target_sd' must be a single positive numeric value.",
        call. = FALSE
      )
    }

    if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau < 0) {
      stop(
        "'tau' must be a single non-negative numeric value.",
        call. = FALSE
      )
    }
  }

  prepared <- prepare_treeangle_data(
    data = data,
    input_type = input_type,
    normalize = normalize_flag,
    alpha = alpha,
    target_sd = target_sd,
    tau = tau
  )

  used_index <- prepared$xy[, 1] > 0 & prepared$xy[, 2] > 0
  n_used <- sum(used_index)

  if (n_used < 2L) {
    stop(
      "Fewer than two points remain after filtering to x > 0 and y > 0. ",
      "Please check the input data or use normalization when generation ",
      "information is available.",
      call. = FALSE
    )
  }

  all_results <- angle_alpha(
    data = prepared$xy,
    alpha = alpha,
    start = c(0, 0),
    return_all = TRUE
  )

  method_table <- data.frame(
    method = names(all_results),
    angle_degrees = vapply(all_results, function(x) unname(x[1, 1]), numeric(1)),
    tan_angle = vapply(all_results, function(x) unname(x[1, 2]), numeric(1)),
    row.names = NULL
  )

  build_treeangle_result(
    result_matrix = all_results[[method]],
    method = method,
    alpha = alpha,
    normalized = prepared$normalized,
    input_type = prepared$input_type,
    n_input = prepared$n_input,
    n_used = n_used,
    n_generations = prepared$n_generations,
    normalization_settings = if (prepared$normalized) {
      list(target_sd = target_sd, tau = tau)
    } else {
      NULL
    },
    method_table = method_table,
    call = call
  )
}

resolve_normalize_flag <- function(normalize, input_type) {
  if (is.null(normalize)) {
    return(input_type %in% c("generation_list", "generation_matrix"))
  }

  if (isTRUE(normalize) && input_type == "xy_matrix") {
    stop(
      "normalize = TRUE requires generation-structured input: either a list of ",
      "generation-specific matrices/data frames or a matrix/data frame with a ",
      "generation column.",
      call. = FALSE
    )
  }

  isTRUE(normalize)
}

prepare_treeangle_data <- function(data,
                                   input_type,
                                   normalize,
                                   alpha,
                                   target_sd,
                                   tau) {
  if (input_type %in% c("generation_list", "generation_matrix")) {
    generation_data <- parse_generation_data(data)

    xy <- if (normalize) {
      do.call(
        rbind,
        normalize_generation_data(
          data_list = generation_data,
          alpha = alpha,
          target_sd = target_sd,
          tau = tau
        )
      )
    } else {
      do.call(rbind, generation_data)
    }

    n_input <- sum(vapply(generation_data, nrow, integer(1)))
    n_generations <- length(generation_data)
  } else {
    xy <- parse_xy_data(data)
    n_input <- nrow(xy)
    n_generations <- NA_integer_
  }

  list(
    xy = xy,
    normalized = normalize,
    input_type = input_type,
    n_input = n_input,
    n_generations = n_generations
  )
}

detect_input_type <- function(data) {
  if (is.list(data) && !is.data.frame(data)) {
    return("generation_list")
  }

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop(
      "'data' must be a list, matrix, data frame, or an object returned by ",
      "'simdata()'.",
      call. = FALSE
    )
  }

  data_matrix <- as.matrix(data)

  if (ncol(data_matrix) < 2L) {
    stop("'data' must contain at least two columns.", call. = FALSE)
  }

  if (ncol(data_matrix) >= 3L) {
    generation_column <- suppressWarnings(as.numeric(data_matrix[, 3]))
    if (all(!is.na(generation_column)) &&
        all(generation_column > 0) &&
        all(generation_column == floor(generation_column))) {
      return("generation_matrix")
    }
  }

  "xy_matrix"
}

ensure_numeric_matrix <- function(x, min_cols, object_name) {
  x_matrix <- as.matrix(x)

  if (ncol(x_matrix) < min_cols) {
    stop(
      object_name,
      " must contain at least ",
      min_cols,
      " columns.",
      call. = FALSE
    )
  }

  x_numeric <- suppressWarnings(
    matrix(
      as.numeric(x_matrix),
      nrow = nrow(x_matrix),
      ncol = ncol(x_matrix),
      dimnames = dimnames(x_matrix)
    )
  )

  if (anyNA(x_numeric[, seq_len(min_cols), drop = FALSE])) {
    stop(
      object_name,
      " contains missing or non-numeric values in the required columns.",
      call. = FALSE
    )
  }

  x_numeric
}

parse_generation_data <- function(data) {
  if (is.list(data) && !is.data.frame(data)) {
    if (length(data) == 0L) {
      stop("'data' is an empty list.", call. = FALSE)
    }

    data_list <- lapply(seq_along(data), function(i) {
      generation_matrix <- ensure_numeric_matrix(
        x = data[[i]],
        min_cols = 2L,
        object_name = paste0("Generation ", i)
      )[, 1:2, drop = FALSE]

      if (nrow(generation_matrix) == 0L) {
        stop("Generation ", i, " contains no observations.", call. = FALSE)
      }

      generation_matrix
    })

    return(data_list)
  }

  data_matrix <- ensure_numeric_matrix(
    x = data,
    min_cols = 3L,
    object_name = "'data'"
  )

  generation_column <- data_matrix[, 3]

  if (any(generation_column < 1) || any(generation_column != floor(generation_column))) {
    stop("The generation column must contain positive integers.", call. = FALSE)
  }

  generations <- sort(unique(generation_column))

  lapply(generations, function(g) {
    generation_matrix <- data_matrix[generation_column == g, 1:2, drop = FALSE]

    if (nrow(generation_matrix) == 0L) {
      stop("Generation ", g, " contains no observations.", call. = FALSE)
    }

    generation_matrix
  })
}

parse_xy_data <- function(data) {
  if (is.list(data) && !is.data.frame(data)) {
    if (length(data) == 0L) {
      stop("'data' is an empty list.", call. = FALSE)
    }

    xy_list <- lapply(seq_along(data), function(i) {
      ensure_numeric_matrix(
        x = data[[i]],
        min_cols = 2L,
        object_name = paste0("Element ", i)
      )[, 1:2, drop = FALSE]
    })

    xy <- do.call(rbind, xy_list)
  } else {
    xy <- ensure_numeric_matrix(
      x = data,
      min_cols = 2L,
      object_name = "'data'"
    )[, 1:2, drop = FALSE]
  }

  if (nrow(xy) < 2L) {
    stop("At least two observations are required.", call. = FALSE)
  }

  xy
}

normalize_generation_data <- function(data_list, alpha, target_sd, tau) {
  lapply(seq_along(data_list), function(i) {
    generation_matrix <- as.matrix(data_list[[i]])
    generation_sd <- sqrt(column_variance_n(generation_matrix))

    if (nrow(generation_matrix) < 2L ||
        any(!is.finite(generation_sd)) ||
        any(generation_sd <= 0)) {
      stop(
        "Each generation must contain at least two non-identical observations ",
        "when normalize = TRUE. In practice, supply generation-wise paired ",
        "observations or increments and do not include a singleton root ",
        "generation.",
        call. = FALSE
      )
    }

    generation_matrix[, 1] <- generation_matrix[, 1] / generation_sd[1] * target_sd
    generation_matrix[, 2] <- generation_matrix[, 2] / generation_sd[2] * target_sd

    harmonic_term <- sum(1 / seq_len(i))
    shift_value <- sqrt(
      -2 * log(1 - alpha) * (target_sd^2) + harmonic_term * tau
    )

    generation_matrix[, 1] <- generation_matrix[, 1] -
      mean(generation_matrix[, 1]) + shift_value
    generation_matrix[, 2] <- generation_matrix[, 2] -
      mean(generation_matrix[, 2]) + shift_value

    generation_matrix
  })
}

column_variance_n <- function(x) {
  x <- as.matrix(x)
  centered <- sweep(x, 2, colMeans(x), FUN = "-")
  colMeans(centered^2)
}

build_treeangle_result <- function(result_matrix,
                                   method,
                                   alpha,
                                   normalized,
                                   input_type,
                                   n_input,
                                   n_used,
                                   n_generations,
                                   normalization_settings,
                                   method_table,
                                   call) {
  rownames(result_matrix) <- c(
    "summary",
    "upper_point",
    "lower_point",
    "upper_line",
    "lower_line"
  )
  colnames(result_matrix) <- c("value_1", "value_2")

  upper_point <- as.numeric(result_matrix["upper_point", ])
  names(upper_point) <- c("x", "y")

  lower_point <- as.numeric(result_matrix["lower_point", ])
  names(lower_point) <- c("x", "y")

  upper_line <- as.numeric(result_matrix["upper_line", ])
  names(upper_line) <- c("slope", "intercept")

  lower_line <- as.numeric(result_matrix["lower_line", ])
  names(lower_line) <- c("slope", "intercept")

  out <- list(
    angle_degrees = unname(result_matrix["summary", 1]),
    tan_angle = unname(result_matrix["summary", 2]),
    upper_point = upper_point,
    lower_point = lower_point,
    upper_line = upper_line,
    lower_line = lower_line,
    method = method,
    alpha = alpha,
    normalized = normalized,
    input_type = input_type,
    n_input = n_input,
    n_used = n_used,
    n_filtered = n_input - n_used,
    n_generations = n_generations,
    normalization_settings = normalization_settings,
    method_table = method_table,
    result_matrix = result_matrix,
    call = call
  )

  class(out) <- "treeangle_result"
  out
}

input_type_label <- function(input_type) {
  switch(
    input_type,
    generation_list = "generation-structured list",
    generation_matrix = "matrix/data frame with generation column",
    xy_matrix = "two-column matrix/data frame",
    input_type
  )
}

#' @export
#' @noRd
print.treeangle_result <- function(x, ...) {
  cat("<treeangle_result>\n")
  cat("  method        : ", x$method, "\n", sep = "")
  cat("  alpha         : ", formatC(x$alpha, digits = 3, format = "f"), "\n", sep = "")
  cat("  normalized    : ", if (x$normalized) "yes" else "no", "\n", sep = "")
  cat("  input_type    : ", input_type_label(x$input_type), "\n", sep = "")

  if (!is.na(x$n_generations)) {
    cat("  generations   : ", x$n_generations, "\n", sep = "")
  }

  cat("  n_input       : ", x$n_input, "\n", sep = "")
  cat("  n_used        : ", x$n_used, "\n", sep = "")
  cat("  n_filtered    : ", x$n_filtered, "\n", sep = "")

  if (x$normalized) {
    cat(
      "  target_sd     : ",
      formatC(x$normalization_settings$target_sd, digits = 4, format = "f"),
      "\n",
      sep = ""
    )
    cat(
      "  tau           : ",
      formatC(x$normalization_settings$tau, digits = 4, format = "f"),
      "\n",
      sep = ""
    )
  }

  cat(
    "  angle_degrees : ",
    formatC(x$angle_degrees, digits = 4, format = "f"),
    "\n",
    sep = ""
  )
  cat(
    "  tan_angle     : ",
    formatC(x$tan_angle, digits = 4, format = "f"),
    "\n",
    sep = ""
  )

  invisible(x)
}

#' @export
#' @noRd
as.matrix.treeangle_result <- function(x, ...) {
  x$result_matrix
}

angle_alpha <- function(data,
                        alpha,
                        start = c(0, 0),
                        return_all = FALSE) {
  data <- as.matrix(data)
  data <- data[(data[, 1] > start[1]) & (data[, 2] > start[2]), , drop = FALSE]

  if (nrow(data) < 2L) {
    stop(
      "Not enough points remain after filtering to x > start[1] and y > start[2].",
      call. = FALSE
    )
  }

  outnum <- floor(nrow(data) * (1 - alpha))

  line_from_point <- function(point) {
    slope <- (point[2] - start[2]) / (point[1] - start[1])
    intercept <- (start[2] * point[1] - point[2] * start[1]) / (point[1] - start[1])
    c(slope, intercept)
  }

  angle_from_point <- function(point) {
    angle <- atan2(point[2] - start[2], point[1] - start[1]) * 180 / pi
    if (angle < 0) angle + 360 else angle
  }

  build_candidate_matrix <- function(angle_value,
                                     upper_point,
                                     lower_point,
                                     upper_line,
                                     lower_line) {
    out <- rbind(
      c(angle_value, tan(angle_value * pi / 180)),
      upper_point,
      lower_point,
      upper_line,
      lower_line
    )
    storage.mode(out) <- "double"
    out
  }

  coef_values <- t(
    vapply(
      seq_len(nrow(data)),
      function(i) line_from_point(data[i, ]),
      numeric(2)
    )
  )

  angle_values <- vapply(
    seq_len(nrow(data)),
    function(i) angle_from_point(data[i, ]),
    numeric(1)
  )

  orderup <- order(-angle_values)
  orderdown <- order(angle_values)

  if (outnum == 0L) {
    upper_index <- orderup[1]
    lower_index <- orderdown[1]
    angle_min <- abs(angle_values[upper_index] - angle_values[lower_index])

    shared_result <- build_candidate_matrix(
      angle_value = angle_min,
      upper_point = data[upper_index, ],
      lower_point = data[lower_index, ],
      upper_line = coef_values[upper_index, ],
      lower_line = coef_values[lower_index, ]
    )

    results <- list(
      mean = shared_result,
      min = shared_result,
      median = shared_result,
      equal = shared_result,
      neigh = shared_result
    )

    return(results)
  }

  candidate_count <- outnum + 1L
  updown <- cbind(
    orderup[seq_len(candidate_count)],
    rev(orderdown[seq_len(candidate_count)])
  )
  anglevs <- abs(angle_values[updown[, 1]] - angle_values[updown[, 2]])

  min_index <- which.min(anglevs)
  min_upper <- updown[min_index, 1]
  min_lower <- updown[min_index, 2]

  angle_min_result <- build_candidate_matrix(
    angle_value = anglevs[min_index],
    upper_point = data[min_upper, ],
    lower_point = data[min_lower, ],
    upper_line = coef_values[min_upper, ],
    lower_line = coef_values[min_lower, ]
  )

  if (candidate_count %% 2L == 1L) {
    middle_index <- (candidate_count + 1L) / 2L
    middle_upper <- updown[middle_index, 1]
    middle_lower <- updown[middle_index, 2]

    angle_equal_result <- build_candidate_matrix(
      angle_value = anglevs[middle_index],
      upper_point = data[middle_upper, ],
      lower_point = data[middle_lower, ],
      upper_line = coef_values[middle_upper, ],
      lower_line = coef_values[middle_lower, ]
    )
  } else {
    left_middle <- candidate_count / 2L
    right_middle <- left_middle + 1L

    upper_point <- colMeans(data[updown[c(left_middle, right_middle), 1], , drop = FALSE])
    lower_point <- colMeans(data[updown[c(left_middle, right_middle), 2], , drop = FALSE])

    angle_equal <- abs(
      mean(angle_values[updown[c(left_middle, right_middle), 1]]) -
        mean(angle_values[updown[c(left_middle, right_middle), 2]])
    )

    angle_equal_result <- build_candidate_matrix(
      angle_value = angle_equal,
      upper_point = upper_point,
      lower_point = lower_point,
      upper_line = line_from_point(upper_point),
      lower_line = line_from_point(lower_point)
    )
  }

  mean_upper_point <- colMeans(data[updown[, 1], , drop = FALSE])
  mean_lower_point <- colMeans(data[updown[, 2], , drop = FALSE])

  angle_mean_result <- build_candidate_matrix(
    angle_value = mean(anglevs),
    upper_point = mean_upper_point,
    lower_point = mean_lower_point,
    upper_line = line_from_point(mean_upper_point),
    lower_line = line_from_point(mean_lower_point)
  )

  median_upper_point <- apply(
    data[updown[, 1], , drop = FALSE],
    2,
    stats::median
  )
  median_lower_point <- apply(
    data[updown[, 2], , drop = FALSE],
    2,
    stats::median
  )

  angle_median_result <- build_candidate_matrix(
    angle_value = stats::median(anglevs),
    upper_point = median_upper_point,
    lower_point = median_lower_point,
    upper_line = line_from_point(median_upper_point),
    lower_line = line_from_point(median_lower_point)
  )

  up_position <- match(min_upper, orderup)
  down_position <- match(min_lower, orderdown)

  neighbor_upper <- orderup[if (up_position == 1L) 1L else up_position - 1L]
  neighbor_lower <- orderdown[if (down_position == 1L) 1L else down_position - 1L]

  neigh_upper_point <- colMeans(data[c(min_upper, neighbor_upper), , drop = FALSE])
  neigh_lower_point <- colMeans(data[c(min_lower, neighbor_lower), , drop = FALSE])

  angle_neigh <- abs(
    mean(angle_values[c(min_upper, neighbor_upper)]) -
      mean(angle_values[c(min_lower, neighbor_lower)])
  )

  angle_neigh_result <- build_candidate_matrix(
    angle_value = angle_neigh,
    upper_point = neigh_upper_point,
    lower_point = neigh_lower_point,
    upper_line = line_from_point(neigh_upper_point),
    lower_line = line_from_point(neigh_lower_point)
  )

  list(
    mean = angle_mean_result,
    min = angle_min_result,
    median = angle_median_result,
    equal = angle_equal_result,
    neigh = angle_neigh_result
  )
}