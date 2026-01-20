#' Calculate angles from 2D data points
#'
#' This function calculates angles from 2D data points using various methods.
#' It filters data based on a starting point and computes angles using the specified method.
#'
#' @param data A matrix or data frame with 2 columns containing x and y coordinates.
#' @param alpha Numeric value between 0 and 1 specifying the quantile threshold (default = 0.95).
#' @param method Character string specifying the angle calculation method. 
#'        Options: "min", "median", "equal", "mean", "neigh" (default = "mean").
#' @param start Numeric vector of length 2 specifying the starting point coordinates (default = c(0,0)).
#'
#' @return A matrix containing angle information, point coordinates, and coefficients for the selected method.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' set.seed(123)
#' sample_data <- matrix(rnorm(200, mean = 5), ncol = 2)
#' 
#' # Calculate angles using default parameters (mean method)
#' result <- calangle(sample_data)
#' 
#' # Calculate angles using min method with custom alpha
#' result_min <- calangle(sample_data, alpha = 0.9, method = "min")
#' }
#'
#' @export
calangle <- function(data, alpha = 0.95, method = "mean", start = c(0,0)) {
  
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Data must be a matrix or data frame")
  }
  
  if (ncol(data) != 2) {
    stop("Data must have exactly 2 columns")
  }
  
  if (alpha <= 0 || alpha >= 1) {
    stop("Alpha must be between 0 and 1")
  }
  
  valid_methods <- c("min", "median", "equal", "mean", "neigh")
  if (!method %in% valid_methods) {
    stop("Method must be one of: 'min', 'median', 'equal', 'mean', 'neigh'")
  }
  
  if (length(start) != 2) {
    stop("Start must be a numeric vector of length 2")
  }
  
  # Convert data to matrix for consistent processing
  data <- as.matrix(data)
  
  # Call the core angle calculation function
  result <- angle_alpha(data, alpha, start)
  
  # Map simplified method names to internal method names
  method_mapping <- c(
    "min" = "anglemin",
    "median" = "anglemedian", 
    "equal" = "angleequal",
    "mean" = "anglemean",
    "neigh" = "angleneigh"
  )
  
  # Return only the selected method
  return(result[[method_mapping[method]]])
}

#' Core angle calculation function
#'
#' This is the internal function that performs the angle calculations.
#'
#' @param data Matrix with 2 columns containing x and y coordinates.
#' @param alpha Numeric value between 0 and 1 for quantile threshold.
#' @param start Numeric vector of length 2 for starting point coordinates.
#'
#' @return A list containing all angle calculation methods.
#'
#' @keywords internal
angle_alpha <- function(data, alpha, start) {
  
  # Filter data based on starting point
  data <- data[((data[,1] > start[1]) & (data[,2] > start[2])), , drop = FALSE]
  
  # Check if enough data points remain
  if (nrow(data) < 2) {
    stop("Not enough data points remaining after filtering by start coordinates.")
  }
  
  # Calculate number of outliers to remove
  outnum <- floor(dim(data)[1] * (1 - alpha))
  
  # Helper function to calculate slope and intercept
  ab <- function(point, start) {
    a <- (point[2] - start[2]) / (point[1] - start[1])
    b <- (start[2] * point[1] - point[2] * start[1]) / (point[1] - start[1])
    return(c(a, b))
  }
  
  # Helper function to calculate angle
  abangle <- function(a, b, point, start) {
    pointtrans <- point - start
    
    if (a >= 0) {
      angle <- atan(a) / pi * 180
    } else if ((a < 0) & (pointtrans[1] < 0)) {
      angle <- atan(a) / pi * 180 + 180
    } else {
      angle <- atan(a) / pi * 180
    }
    return(angle)
  }
  
  # Initialize vectors
  Tan <- c()
  Angle <- c()
  Coef <- c()
  
  # Loop through data points
  for (i in 1:dim(data)[1]) {
    point <- data[i, ]
    coef <- ab(point, start)
    Tan <- c(Tan, coef[1])
    Coef <- rbind(Coef, coef)
    
    angle <- abangle(coef[1], coef[2], point, start)
    Angle <- c(Angle, angle)
  }
  
  # Helper function for comparisons
  comple <- function(outnum, Angle, orderup, orderdown) {
    updown <- c()
    anglevs <- c()
    for (i in 1:(outnum + 1)) {
      up <- orderup[i]
      down <- orderdown[outnum + 2 - i]
      updown <- rbind(updown, c(up, down))
      anglevs <- c(anglevs, abs(Angle[up] - Angle[down]))
    }
    
    return(list(updown, anglevs))
  }
  
  orderup <- order(-Angle)
  orderdown <- order(Angle)
  
  # Core logic handling different cases for outnum
  if (outnum == 0) {
    whereup <- orderup[1]
    wheredown <- orderdown[1]
    anglemin <- abs(Angle[whereup] - Angle[wheredown])
    tanup <- Tan[whereup]
    tandown <- Tan[wheredown]
    tan2 <- (tanup - tandown) / (1 + tanup * tandown)
    pointup <- data[whereup, ]
    coefup <- Coef[whereup, ]
    pointdown <- data[wheredown, ]
    coefdown <- Coef[wheredown, ]
    
    angleminlist <- rbind(c(anglemin, tan2), pointup, pointdown, coefup, coefdown)
    angleequallist <- angleminlist
    anglemeanlist <- angleminlist
    angleneighlist <- angleminlist
    anglemedianlist <- angleminlist
  } else {
    updownanglevs <- comple(outnum, Angle, orderup, orderdown)
    updown <- updownanglevs[[1]]
    anglevs <- updownanglevs[[2]]
    
    where <- order(anglevs)[1]
    whereup <- updown[where, 1]
    wheredown <- updown[where, 2]
    anglemin <- anglevs[where]
    tanup <- Tan[whereup]
    tandown <- Tan[wheredown]
    tan2 <- (tanup - tandown) / (1 + tanup * tandown)
    pointup <- data[whereup, ]
    coefup <- Coef[whereup, ]
    pointdown <- data[wheredown, ]
    coefdown <- Coef[wheredown, ]
    angleminlist <- rbind(c(anglemin, tan2), pointup, pointdown, coefup, coefdown)
    
    if (outnum == 1) {
      angleequallist <- angleminlist
    } else if (outnum %% 2 == 0) {
      where <- ceiling(outnum / 2)
      whereup <- updown[where, 1]
      wheredown <- updown[where, 2]
      angleequal <- anglevs[where]
      tanup <- Tan[whereup]
      tandown <- Tan[wheredown]
      tan2 <- (tanup - tandown) / (1 + tanup * tandown)
      pointup <- data[whereup, ]
      coefup <- Coef[whereup, ]
      pointdown <- data[wheredown, ]
      coefdown <- Coef[wheredown, ]
      angleequallist <- rbind(c(angleequal, tan2), pointup, pointdown, coefup, coefdown)
    } else {
      where <- floor(outnum / 2)
      whereup <- updown[where, 1]
      whereupup <- updown[where + 1, 1]
      wheredown <- updown[where, 2]
      wheredowndown <- updown[where + 1, 2]
      pointup <- apply(data[c(whereup, whereupup), , drop = FALSE], 2, mean)
      coefup <- ab(pointup, start)
      pointdown <- apply(data[c(wheredown, wheredowndown), , drop = FALSE], 2, mean)
      coefdown <- ab(pointdown, start)
      angleequal <- abs(mean(Angle[c(whereup, whereupup)]) - mean(Angle[c(wheredown, wheredowndown)]))
      tan2 <- tan(angleequal)
      angleequallist <- rbind(c(angleequal, tan2), pointup, pointdown, coefup, coefdown)
    }
    
    anglemean <- mean(anglevs)
    tan2 <- tan(anglemean)
    pointup <- apply(data[updown[, 1], , drop = FALSE], 2, mean)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[updown[, 2], , drop = FALSE], 2, mean)
    coefdown <- ab(pointdown, start)
    anglemeanlist <- rbind(c(anglemean, tan2), pointup, pointdown, coefup, coefdown)
    
    anglemedian <- median(anglevs)
    tan2 <- tan(anglemedian)
    pointup <- apply(data[updown[, 1], , drop = FALSE], 2, median)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[updown[, 2], , drop = FALSE], 2, median)
    coefdown <- ab(pointdown, start)
    anglemedianlist <- rbind(c(anglemedian, tan2), pointup, pointdown, coefup, coefdown)
    
    where <- order(anglevs)[1]
    whereup <- updown[where, 1]
    upnum <- which(orderup == whereup)
    whereupup <- orderup[ifelse(upnum == 1, upnum, upnum - 1)]
    wheredown <- updown[where, 2]
    downnum <- which(orderdown == wheredown)
    wheredowndown <- orderdown[ifelse(downnum == 1, downnum, downnum - 1)]
    pointup <- apply(data[c(whereup, whereupup), , drop = FALSE], 2, mean)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[c(wheredown, wheredowndown), , drop = FALSE], 2, mean)
    coefdown <- ab(pointdown, start)
    angleneigh <- abs(mean(Angle[c(whereup, whereupup)]) - mean(Angle[c(wheredown, wheredowndown)]))
    tan2 <- tan(angleneigh)
    angleneighlist <- rbind(c(angleneigh, tan2), pointup, pointdown, coefup, coefdown)
  }
  
  out <- list(anglemin = angleminlist, 
              anglemedian = anglemedianlist, 
              angleequal = angleequallist, 
              anglemean = anglemeanlist, 
              angleneigh = angleneighlist)
  
  return(out)
}