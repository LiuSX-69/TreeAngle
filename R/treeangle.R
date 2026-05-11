#' Compute the Treeangle statistic
#'
#' `treeangle()` computes the sample Treeangle statistic from paired
#' two-dimensional observations or increments.
#'
#' For generation-structured input, the normalization step in the manuscript
#' can be applied before calculating the angle. The default behavior follows the
#' TD Delta theta algorithm: generation-structured data are normalized
#' automatically, while direct two-column data are not normalized.
#'
#' The output is deliberately kept close to the original script:
#'
#' - `method = "mean"` (default) returns the average of all valid candidate angles;
#' - `method = "min"` returns the smallest valid candidate angle;
#' - `method = "median"` returns the median of all valid candidate angles;
#' - `method = "equal"` returns the candidate chosen by the equal-allocation rule from the original implementation;
#' - `method = "neigh"` returns the candidate chosen by the neighborhood-adjusted rule from the original implementation;
#' - `method = "all"` returns the full list
#'   `min`, `median`, `equal`, `mean`, `neigh`.
#'
#' @param data Input data. Supported formats:
#' \itemize{
#'   \item a data frame or matrix with columns `generation`, `x`, `y`;
#'   \item a data frame or matrix with columns `x`, `y`, `generation`;
#'   \item a two-column matrix/data frame of paired observations or increments;
#'   \item a list of generation-specific two-column matrices/data frames;
#'   \item an old-style list `list(treeresidual1, treeresidual2)`, where each
#'   component is a generation list.
#' }
#' @param alpha Numeric scalar in `(0, 1)`. Coverage proportion. For example,
#' `alpha = 0.95` means that 95 percent of points are kept inside the angle.
#' @param method Character string. One of `"mean"`, `"min"`, `"median"`,
#' `"equal"`, `"neigh"`, or `"all"`. The original names `"anglemean"`,
#' `"anglemin"`, `"anglemedian"`, `"angleequal"`, and `"angleneigh"` are also
#' accepted.
#' @param normalize Logical scalar or `NULL`. If `NULL`, normalization is applied
#' automatically to generation-structured input and skipped for direct
#' two-column input.
#' @param target_sd Positive numeric scalar. Target marginal standard deviation
#' used in normalization. The manuscript simulation uses `target_sd = 1`.
#' To reproduce the old script target variance `0.0001`, use
#' `target_sd = 0.01`.
#' @param tau Non-negative numeric scalar. Generation shift parameter in the
#' normalization step. Default is `0.1`.
#' @param start Numeric vector of length two. Fixed intersect point of the two
#' rays. Default is `c(0, 0)`.
#'
#' @return A plain matrix or list.
#'
#' If `method != "all"`, the return value is a 5 x 2 numeric matrix.
#' The row names are:
#' \itemize{
#'   \item `angle`: angle summary;
#'   \item `pointup`: point on the upper boundary ray;
#'   \item `pointdown`: point on the lower boundary ray;
#'   \item `coefup`: slope and intercept of the upper boundary line;
#'   \item `coefdown`: slope and intercept of the lower boundary line.
#' }
#'
#' The column names are `deg_x_slope` and `tan_y_intercept`.
#' Their meanings depend on the row:
#' \itemize{
#'   \item row `angle`: angle in degrees and tangent of the angle;
#'   \item rows `pointup` and `pointdown`: x-coordinate and y-coordinate;
#'   \item rows `coefup` and `coefdown`: slope and intercept.
#' }
#'
#' If `method = "all"`, the return value is a list with names
#' `anglemin`, `anglemedian`, `angleequal`, `anglemean`, and `angleneigh`.
#'
#' @examples
#' dat <- simdata(seed = 1)
#'
#' ## Paper default: normalized TD Delta theta, using mean candidate angle
#' treeangle(dat, alpha = 0.95, method = "mean")
#'
#' ## Return all five original summaries
#' treeangle(dat, alpha = 0.95, method = "all")
#'
#' ## Without normalization, corresponding to the preliminary angle calculation
#' treeangle(dat, alpha = 0.95, method = "mean", normalize = FALSE)
#'
#' @export
treeangle <- function(data,
                      alpha = 0.95,
                      method = "mean",
                      normalize = NULL,
                      target_sd = 1,
                      tau = 0.1,
                      start = c(0, 0)) {
  validate_treeangle_alpha(alpha)

  method <- match_treeangle_method(method)

  if (
    !is.null(normalize) &&
      !(is.logical(normalize) &&
        length(normalize) == 1L &&
        !is.na(normalize))
  ) {
    stop("'normalize' must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  if (
    !is.numeric(start) ||
      length(start) != 2L ||
      anyNA(start) ||
      any(!is.finite(start))
  ) {
    stop("'start' must be a finite numeric vector of length 2.", call. = FALSE)
  }

  start <- as.numeric(start)

  if (inherits(data, "treeangle_simulation")) {
    data <- data$data
  }

  input_type <- detect_input_type(data)

  normalize_flag <- resolve_normalize_flag(
    normalize = normalize,
    input_type = input_type
  )

  if (normalize_flag) {
    if (
      !is.numeric(target_sd) ||
        length(target_sd) != 1L ||
        is.na(target_sd) ||
        !is.finite(target_sd) ||
        target_sd <= 0
    ) {
      stop("'target_sd' must be a single positive numeric value.", call. = FALSE)
    }

    if (
      !is.numeric(tau) ||
        length(tau) != 1L ||
        is.na(tau) ||
        !is.finite(tau) ||
        tau < 0
    ) {
      stop("'tau' must be a single non-negative numeric value.", call. = FALSE)
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

  out <- anglealpha(
    data = prepared$xy,
    alpha = alpha,
    start = start
  )

  if (identical(method, "all")) {
    return(out)
  }

  out[[method]]
}

validate_treeangle_alpha <- function(alpha) {
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop(
      "'alpha' must be a single numeric value strictly between 0 and 1.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

match_treeangle_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("'method' must be a single character string.", call. = FALSE)
  }

  method_map <- c(
    min = "anglemin",
    median = "anglemedian",
    equal = "angleequal",
    mean = "anglemean",
    neigh = "angleneigh",
    anglemin = "anglemin",
    anglemedian = "anglemedian",
    angleequal = "angleequal",
    anglemean = "anglemean",
    angleneigh = "angleneigh",
    all = "all"
  )

  if (!(method %in% names(method_map))) {
    stop(
      "'method' must be one of: mean, min, median, equal, neigh, all.",
      call. = FALSE
    )
  }

  unname(method_map[[method]])
}

resolve_normalize_flag <- function(normalize, input_type) {
  generation_types <- c(
    "generation_list",
    "generation_matrix",
    "paired_generation_list"
  )

  if (is.null(normalize)) {
    return(input_type %in% generation_types)
  }

  if (isTRUE(normalize) && input_type == "xy_matrix") {
    stop(
      "normalize = TRUE requires generation-structured input.",
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
  generation_types <- c(
    "generation_list",
    "generation_matrix",
    "paired_generation_list"
  )

  if (input_type %in% generation_types) {
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
  } else {
    xy <- parse_xy_data(data)
  }

  list(xy = xy)
}

has_named_columns <- function(x, cols) {
  nms <- colnames(x)
  !is.null(nms) && all(cols %in% nms)
}

select_numeric_columns <- function(x, cols, object_name) {
  x_matrix <- as.matrix(x)

  if (is.character(cols)) {
    missing_cols <- setdiff(cols, colnames(x_matrix))

    if (length(missing_cols) > 0L) {
      stop(
        object_name,
        " is missing required column(s): ",
        paste(missing_cols, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    selected <- x_matrix[, cols, drop = FALSE]
  } else {
    if (ncol(x_matrix) < max(cols)) {
      stop(
        object_name,
        " must contain at least ",
        max(cols),
        " columns.",
        call. = FALSE
      )
    }

    selected <- x_matrix[, cols, drop = FALSE]
  }

  out <- suppressWarnings(
    matrix(
      as.numeric(selected),
      nrow = nrow(selected),
      ncol = ncol(selected)
    )
  )

  if (anyNA(out) || any(!is.finite(out))) {
    stop(
      object_name,
      " contains missing, infinite, or non-numeric values in the required ",
      "columns.",
      call. = FALSE
    )
  }

  storage.mode(out) <- "double"
  out
}

is_paired_generation_list <- function(data) {
  is.list(data) &&
    !is.data.frame(data) &&
    length(data) == 2L &&
    is.list(data[[1]]) &&
    is.list(data[[2]]) &&
    !is.data.frame(data[[1]]) &&
    !is.data.frame(data[[2]]) &&
    length(data[[1]]) > 0L &&
    length(data[[1]]) == length(data[[2]])
}

detect_input_type <- function(data) {
  if (is_paired_generation_list(data)) {
    return("paired_generation_list")
  }

  if (is.list(data) && !is.data.frame(data)) {
    return("generation_list")
  }

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop(
      "'data' must be a list, matrix, data frame, or an object returned by ",
      "simdata().",
      call. = FALSE
    )
  }

  data_matrix <- as.matrix(data)

  if (ncol(data_matrix) < 2L) {
    stop("'data' must contain at least two columns.", call. = FALSE)
  }

  if (has_named_columns(data, c("x", "y", "generation"))) {
    return("generation_matrix")
  }

  if (has_named_columns(data, c("x", "y"))) {
    return("xy_matrix")
  }

  if (ncol(data_matrix) >= 3L) {
    generation_column <- suppressWarnings(as.numeric(data_matrix[, 3]))

    if (
      all(!is.na(generation_column)) &&
        all(is.finite(generation_column)) &&
        all(generation_column > 0) &&
        all(generation_column == floor(generation_column))
    ) {
      return("generation_matrix")
    }
  }

  "xy_matrix"
}

parse_generation_data <- function(data) {
  if (is_paired_generation_list(data)) {
    return(parse_paired_generation_list(data))
  }

  if (is.list(data) && !is.data.frame(data)) {
    if (length(data) == 0L) {
      stop("'data' is an empty list.", call. = FALSE)
    }

    data_list <- lapply(seq_along(data), function(i) {
      cols <- if (has_named_columns(data[[i]], c("x", "y"))) {
        c("x", "y")
      } else {
        1:2
      }

      generation_matrix <- select_numeric_columns(
        x = data[[i]],
        cols = cols,
        object_name = paste0("Generation ", i)
      )[, 1:2, drop = FALSE]

      colnames(generation_matrix) <- c("x", "y")

      if (nrow(generation_matrix) == 0L) {
        stop("Generation ", i, " contains no observations.", call. = FALSE)
      }

      generation_matrix
    })

    names(data_list) <- paste0("generation_", seq_along(data_list))
    return(data_list)
  }

  cols <- if (has_named_columns(data, c("x", "y", "generation"))) {
    c("x", "y", "generation")
  } else {
    1:3
  }

  data_matrix <- select_numeric_columns(
    x = data,
    cols = cols,
    object_name = "'data'"
  )

  colnames(data_matrix) <- c("x", "y", "generation")

  generation_column <- data_matrix[, "generation"]

  if (
    any(generation_column < 1) ||
      any(generation_column != floor(generation_column))
  ) {
    stop("The generation column must contain positive integers.", call. = FALSE)
  }

  generations <- sort(unique(as.integer(generation_column)))
  expected_generations <- seq_len(max(generations))

  if (!identical(generations, expected_generations)) {
    missing_generations <- setdiff(expected_generations, generations)

    stop(
      "The generation column must contain consecutive positive integers ",
      "starting from 1. Missing generation(s): ",
      paste(missing_generations, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  out <- lapply(generations, function(g) {
    generation_matrix <- data_matrix[
      generation_column == g,
      c("x", "y"),
      drop = FALSE
    ]

    if (nrow(generation_matrix) == 0L) {
      stop("Generation ", g, " contains no observations.", call. = FALSE)
    }

    generation_matrix
  })

  names(out) <- paste0("generation_", generations)
  out
}

parse_paired_generation_list <- function(data) {
  x_list <- data[[1]]
  y_list <- data[[2]]

  out <- lapply(seq_along(x_list), function(i) {
    x_i <- suppressWarnings(as.numeric(x_list[[i]]))
    y_i <- suppressWarnings(as.numeric(y_list[[i]]))

    if (length(x_i) != length(y_i)) {
      stop(
        "Generation ",
        i,
        " has unequal x and y lengths.",
        call. = FALSE
      )
    }

    if (length(x_i) == 0L) {
      stop("Generation ", i, " contains no observations.", call. = FALSE)
    }

    if (
      anyNA(x_i) ||
        anyNA(y_i) ||
        any(!is.finite(x_i)) ||
        any(!is.finite(y_i))
    ) {
      stop(
        "Generation ",
        i,
        " contains missing, infinite, or non-numeric values.",
        call. = FALSE
      )
    }

    cbind(x = x_i, y = y_i)
  })

  names(out) <- paste0("generation_", seq_along(out))
  out
}

parse_xy_data <- function(data) {
  cols <- if (has_named_columns(data, c("x", "y"))) {
    c("x", "y")
  } else {
    1:2
  }

  xy <- select_numeric_columns(
    x = data,
    cols = cols,
    object_name = "'data'"
  )[, 1:2, drop = FALSE]

  colnames(xy) <- c("x", "y")

  if (nrow(xy) < 2L) {
    stop("At least two observations are required.", call. = FALSE)
  }

  xy
}

normalize_generation_data <- function(data_list, alpha, target_sd, tau) {
  if (length(data_list) == 0L) {
    stop("'data_list' must contain at least one generation.", call. = FALSE)
  }

  normalized_list <- lapply(seq_along(data_list), function(i) {
    generation_matrix <- as.matrix(data_list[[i]])

    if (nrow(generation_matrix) < 2L) {
      stop(
        "Generation ",
        i,
        " contains fewer than two observations. ",
        "Per-generation normalization requires at least two observations.",
        call. = FALSE
      )
    }

    generation_sd <- sqrt(column_variance_n(generation_matrix))

    if (any(!is.finite(generation_sd)) || any(generation_sd <= 0)) {
      stop(
        "Generation ",
        i,
        " has zero or invalid empirical marginal standard deviation.",
        call. = FALSE
      )
    }

    generation_matrix[, 1] <- generation_matrix[, 1] /
      generation_sd[1] *
      target_sd
    generation_matrix[, 2] <- generation_matrix[, 2] /
      generation_sd[2] *
      target_sd

    harmonic_term <- sum(1 / seq_len(i))

    shift_value <- sqrt(
      -2 * log(1 - alpha) * target_sd^2 + harmonic_term * tau
    )

    generation_matrix[, 1] <- generation_matrix[, 1] -
      mean(generation_matrix[, 1]) +
      shift_value

    generation_matrix[, 2] <- generation_matrix[, 2] -
      mean(generation_matrix[, 2]) +
      shift_value

    colnames(generation_matrix) <- c("x", "y")
    generation_matrix
  })

  names(normalized_list) <- names(data_list)
  normalized_list
}

column_variance_n <- function(x) {
  x <- as.matrix(x)
  centered <- sweep(x, 2, colMeans(x), FUN = "-")
  colMeans(centered^2)
}

anglealpha <- function(data,
                       alpha = 0.95,
                       start = c(0, 0)) {
  validate_treeangle_alpha(alpha)

  if (
    !is.numeric(start) ||
      length(start) != 2L ||
      anyNA(start) ||
      any(!is.finite(start))
  ) {
    stop("'start' must be a finite numeric vector of length 2.", call. = FALSE)
  }

  start <- as.numeric(start)

  data <- as.matrix(data)

  if (ncol(data) < 2L) {
    stop("'data' must contain at least two columns.", call. = FALSE)
  }

  data <- data[, 1:2, drop = FALSE]

  data <- suppressWarnings(
    matrix(
      as.numeric(data),
      nrow = nrow(data),
      ncol = 2
    )
  )

  if (anyNA(data) || any(!is.finite(data))) {
    stop(
      "'data' contains missing, infinite, or non-numeric values.",
      call. = FALSE
    )
  }

  data <- data[
    (data[, 1] > start[1]) & (data[, 2] > start[2]),
    ,
    drop = FALSE
  ]

  if (nrow(data) < 2L) {
    stop(
      "Fewer than two points remain after filtering by start.",
      call. = FALSE
    )
  }

  out_num <- floor(nrow(data) * (1 - alpha))

  ab <- function(point, start) {
    a <- (point[2] - start[2]) / (point[1] - start[1])
    b <- (start[2] * point[1] - point[2] * start[1]) /
      (point[1] - start[1])
    c(a, b)
  }

  abangle <- function(a, b, point, start) {
    pointtrans <- point - start

    if (a >= 0) {
      angle <- atan(a) / pi * 180
    } else if ((a < 0) && (pointtrans[1] < 0)) {
      angle <- atan(a) / pi * 180 + 180
    } else {
      angle <- atan(a) / pi * 180
    }

    angle
  }

  tan_from_degrees <- function(angle_degrees) {
    tan(angle_degrees * pi / 180)
  }

  tan_between_slopes <- function(tanup, tandown) {
    (tanup - tandown) / (1 + tanup * tandown)
  }

  legacy_matrix <- function(angle_value,
                            tan_value,
                            pointup,
                            pointdown,
                            coefup,
                            coefdown) {
    out <- matrix(
      as.numeric(c(
        angle_value, tan_value,
        pointup[1], pointup[2],
        pointdown[1], pointdown[2],
        coefup[1], coefup[2],
        coefdown[1], coefdown[2]
      )),
      nrow = 5,
      ncol = 2,
      byrow = TRUE,
      dimnames = list(
        c(
          "angle",
          "pointup",
          "pointdown",
          "coefup",
          "coefdown"
        ),
        c(
          "deg_x_slope",
          "tan_y_intercept"
        )
      )
    )

    storage.mode(out) <- "double"
    out
  }

  coef_list <- lapply(seq_len(nrow(data)), function(i) ab(data[i, ], start))
  coef_mat <- do.call(rbind, coef_list)
  tan_vec <- coef_mat[, 1]
  angle <- vapply(
    seq_len(nrow(data)),
    function(i) abangle(coef_mat[i, 1], coef_mat[i, 2], data[i, ], start),
    numeric(1)
  )

  compute_candidate_pairs <- function(out_num, angle, order_up, order_down) {
    updown <- NULL
    anglevs <- numeric(0)

    for (i in seq_len(out_num + 1L)) {
      up <- order_up[i]
      down <- order_down[out_num + 2L - i]

      updown <- rbind(updown, c(up, down))
      anglevs <- c(anglevs, abs(angle[up] - angle[down]))
    }

    list(updown, anglevs)
  }

  order_up <- order(-angle)
  order_down <- order(angle)

  if (out_num == 0L) {
    whereup <- order_up[1]
    wheredown <- order_down[1]

    anglemin <- abs(angle[whereup] - angle[wheredown])

    tanup <- tan_vec[whereup]
    tandown <- tan_vec[wheredown]
    tan2 <- tan_between_slopes(tanup, tandown)

    pointup <- data[whereup, ]
    coefup <- coef_mat[whereup, ]

    pointdown <- data[wheredown, ]
    coefdown <- coef_mat[wheredown, ]

    angleminlist <- legacy_matrix(
      angle_value = anglemin,
      tan_value = tan2,
      pointup = pointup,
      pointdown = pointdown,
      coefup = coefup,
      coefdown = coefdown
    )

    angleequallist <- angleminlist
    anglemeanlist <- angleminlist
    angleneighlist <- angleminlist
    anglemedianlist <- angleminlist
  } else {
    updownanglevs <- compute_candidate_pairs(out_num, angle, order_up, order_down)
    updown <- updownanglevs[[1]]
    anglevs <- updownanglevs[[2]]

    where <- order(anglevs)[1]
    whereup <- updown[where, 1]
    wheredown <- updown[where, 2]

    anglemin <- anglevs[where]

    tanup <- tan_vec[whereup]
    tandown <- tan_vec[wheredown]
    tan2 <- tan_between_slopes(tanup, tandown)

    pointup <- data[whereup, ]
    coefup <- coef_mat[whereup, ]

    pointdown <- data[wheredown, ]
    coefdown <- coef_mat[wheredown, ]

    angleminlist <- legacy_matrix(
      angle_value = anglemin,
      tan_value = tan2,
      pointup = pointup,
      pointdown = pointdown,
      coefup = coefup,
      coefdown = coefdown
    )

    if (out_num == 1L) {
      angleequallist <- angleminlist
    } else if (out_num %% 2L == 0L) {
      where <- ceiling(out_num / 2)
      whereup <- updown[where, 1]
      wheredown <- updown[where, 2]

      angleequal <- anglevs[where]

      tanup <- tan_vec[whereup]
      tandown <- tan_vec[wheredown]
      tan2 <- tan_between_slopes(tanup, tandown)

      pointup <- data[whereup, ]
      coefup <- coef_mat[whereup, ]

      pointdown <- data[wheredown, ]
      coefdown <- coef_mat[wheredown, ]

      angleequallist <- legacy_matrix(
        angle_value = angleequal,
        tan_value = tan2,
        pointup = pointup,
        pointdown = pointdown,
        coefup = coefup,
        coefdown = coefdown
      )
    } else {
      where <- floor(out_num / 2)

      whereup <- updown[where, 1]
      whereupup <- updown[where + 1L, 1]

      wheredown <- updown[where, 2]
      wheredowndown <- updown[where + 1L, 2]

      pointup <- apply(data[c(whereup, whereupup), , drop = FALSE], 2, mean)
      coefup <- ab(pointup, start)

      pointdown <- apply(
        data[c(wheredown, wheredowndown), , drop = FALSE],
        2,
        mean
      )
      coefdown <- ab(pointdown, start)

      angleequal <- abs(
        mean(angle[c(whereup, whereupup)]) -
          mean(angle[c(wheredown, wheredowndown)])
      )

      tan2 <- tan_from_degrees(angleequal)

      angleequallist <- legacy_matrix(
        angle_value = angleequal,
        tan_value = tan2,
        pointup = pointup,
        pointdown = pointdown,
        coefup = coefup,
        coefdown = coefdown
      )
    }

    anglemean <- mean(anglevs)
    tan2 <- tan_from_degrees(anglemean)

    pointup <- apply(data[updown[, 1], , drop = FALSE], 2, mean)
    coefup <- ab(pointup, start)

    pointdown <- apply(data[updown[, 2], , drop = FALSE], 2, mean)
    coefdown <- ab(pointdown, start)

    anglemeanlist <- legacy_matrix(
      angle_value = anglemean,
      tan_value = tan2,
      pointup = pointup,
      pointdown = pointdown,
      coefup = coefup,
      coefdown = coefdown
    )

    anglemedian <- stats::median(anglevs)
    tan2 <- tan_from_degrees(anglemedian)

    pointup <- apply(data[updown[, 1], , drop = FALSE], 2, stats::median)
    coefup <- ab(pointup, start)

    pointdown <- apply(data[updown[, 2], , drop = FALSE], 2, stats::median)
    coefdown <- ab(pointdown, start)

    anglemedianlist <- legacy_matrix(
      angle_value = anglemedian,
      tan_value = tan2,
      pointup = pointup,
      pointdown = pointdown,
      coefup = coefup,
      coefdown = coefdown
    )

    where <- order(anglevs)[1]

    whereup <- updown[where, 1]
    upnum <- which(order_up == whereup)
    whereupup <- order_up[ifelse(upnum == 1L, upnum, upnum - 1L)]

    wheredown <- updown[where, 2]
    downnum <- which(order_down == wheredown)
    wheredowndown <- order_down[ifelse(downnum == 1L, downnum, downnum - 1L)]

    pointup <- apply(data[c(whereup, whereupup), , drop = FALSE], 2, mean)
    coefup <- ab(pointup, start)

    pointdown <- apply(
      data[c(wheredown, wheredowndown), , drop = FALSE],
      2,
      mean
    )
    coefdown <- ab(pointdown, start)

    angleneigh <- abs(
      mean(angle[c(whereup, whereupup)]) -
        mean(angle[c(wheredown, wheredowndown)])
    )

    tan2 <- tan_from_degrees(angleneigh)

    angleneighlist <- legacy_matrix(
      angle_value = angleneigh,
      tan_value = tan2,
      pointup = pointup,
      pointdown = pointdown,
      coefup = coefup,
      coefdown = coefdown
    )
  }

  list(
    anglemin = angleminlist,
    anglemedian = anglemedianlist,
    angleequal = angleequallist,
    anglemean = anglemeanlist,
    angleneigh = angleneighlist
  )
}