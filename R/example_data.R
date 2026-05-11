#' Read the example TreeAngle data
#'
#' This function reads a small example CSV file shipped with the package.
#' The file was generated using [simdata()] with a fixed random seed.
#'
#' The returned data frame shows the recommended real-data input format:
#' one row per paired observation or increment, with columns `generation`,
#' `x`, and `y`.
#'
#' @return A data frame with columns `generation`, `x`, and `y`.
#'
#' @examples
#' dat <- treeangle_example_data()
#' head(dat)
#' table(dat$generation)
#'
#' res <- treeangle(
#'   dat,
#'   alpha = 0.95,
#'   method = "mean",
#'   normalize = TRUE
#' )
#'
#' print(res)
#'
#' @export
treeangle_example_data <- function() {
  path <- system.file(
    "extdata",
    "treeangle_example.csv",
    package = "TreeAngle"
  )

  if (!nzchar(path)) {
    stop(
      "Example data file was not found. Please reinstall TreeAngle.",
      call. = FALSE
    )
  }

  utils::read.csv(path, stringsAsFactors = FALSE)
}