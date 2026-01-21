#' Calculate angles from Tree-Structured Data with Generation-Dependent Normalization
#'
#' This function calculates the Minimal Enclosing Angle from 2D data points.
#' It includes a specialized normalization step designed for tree-structured data,
#' strictly following the TD-DeltaTheta algorithm (Algorithm 1 in the paper).
#' This involves scaling based on the last generation's variance and shifting means
#' dynamically based on the generation depth.
#'
#' @param data Input data. Can be one of two formats:
#'        \itemize{
#'          \item **List**: A list of length \eqn{T}, where \code{data[[i]]} is a matrix or data frame
#'                of 2 columns (X, Y) representing the data for generation \eqn{i}.
#'                (This matches the structure of the simulation code in the original study).
#'          \item **Matrix/Data Frame**: Must have 3 columns (X, Y, Generation) if \code{normalize = TRUE}.
#'                The 3rd column must contain integers indicating the generation index.
#'        }
#' @param alpha Numeric value between 0 and 1 specifying the quantile threshold (default = 0.95).
#'              This determines the size of the quantile ellipse (filtering top 5\% outliers).
#' @param method Character string specifying the angle calculation method.
#'        Options: "min", "median", "equal", "mean", "neigh" (default = "mean").
#'        The paper recommends "mean" for robust estimation.
#' @param start Numeric vector of length 2 specifying the starting point coordinates (default = c(0,0)).
#' @param normalize Logical. Whether to perform the generation-dependent normalization (default = TRUE).
#'        If TRUE, the function applies scaling and shifting logic to standardize distributions across generations.
#' @param target_sd Numeric value. The target standard deviation after normalization (default = 0.01).
#'        Corresponds to the \code{sqrt(0.0001)} term in the original algorithm.
#' @param tau Numeric value. A hyperparameter controlling the step size of the mean shift across generations (default = 0.1).
#'        Larger values separate the generations further apart from the origin.
#'        Corresponds to the coefficient for the harmonic series sum \eqn{\sum (1/k)} in the algorithm.
#'
#' @return A matrix containing angle information, point coordinates, and coefficients for the selected method.
#'
#' @examples
#' \dontrun{
#' # --- Scenario 1: Input as a List (Matches original simulation) ---
#' set.seed(123)
#' data_list <- list()
#' # Gen 1
#' data_list[[1]] <- matrix(rnorm(100), ncol=2)
#' # Gen 2
#' data_list[[2]] <- matrix(rnorm(100, mean=2), ncol=2)
#'
#' # Calculate with normalization (Auto-detects generation by list index)
#' res_list <- treeangle(data_list, normalize = TRUE, tau = 0.1)
#'
#' # --- Scenario 2: Input as a Matrix with Generation Column ---
#' data_mat <- rbind(
#'   cbind(data_list[[1]], 1), # Col 3 is Gen
#'   cbind(data_list[[2]], 2)
#' )
#'
#' res_mat <- treeangle(data_mat, normalize = TRUE)
#' }
#'
#' @export
treeangle <- function(data, alpha = 0.95, method = "mean", start = c(0,0),
                     normalize = TRUE, target_sd = 0.01, tau = 0.1) {
    
  data_list <- list()
  
  if (is.list(data) && !is.data.frame(data)) {
    
    data_list <- lapply(data, as.matrix)
    
  } else {

    if (!is.matrix(data) && !is.data.frame(data)) {
      stop("Data must be a matrix, data frame, or a list of matrices.")
    }
    data <- as.matrix(data)
    
    if (normalize) {
      if (ncol(data) < 3) {
        stop("When 'normalize = TRUE' and input is a matrix, it must have at least 3 columns: (X, Y, Generation).")
      }
      gen_col <- data[, 3]
      gens <- sort(unique(gen_col))
      if (min(gens) < 1) stop("Generation column must contain positive integers.")
      
      data_list <- lapply(gens, function(g) {
        data[gen_col == g, 1:2, drop = FALSE]
      })
    } else {
      
      if (ncol(data) < 2) stop("Data must have at least 2 columns (X, Y).")
      data_list <- list(data[, 1:2, drop = FALSE])
    }
  }
  
  # ---  Normalization  ---
  if (normalize) {
    T_max <- length(data_list)
    if (T_max < 1) stop("Data is empty.")
    
    last_gen_data <- data_list[[T_max]]
    if (nrow(last_gen_data) < 2) stop("Last generation has insufficient data to calculate variance.")
    
    sigma_ref_x <- sd(last_gen_data[,1])
    sigma_ref_y <- sd(last_gen_data[,2])
    
    # Safety check for constant data
    if(sigma_ref_x == 0) sigma_ref_x <- 1
    if(sigma_ref_y == 0) sigma_ref_y <- 1
    
    # Step B: Loop through generations i = 1 to T
    for (i in 1:T_max) {
      gen_data <- data_list[[i]]
      
      # 1. Centering (Subtract current generation mean)
      mu_x <- mean(gen_data[,1])
      mu_y <- mean(gen_data[,2])
      
      scaled_x <- (gen_data[,1] - mu_x) / sigma_ref_x * target_sd
      scaled_y <- (gen_data[,2] - mu_y) / sigma_ref_y * target_sd

      term1 <- -2 * log(1 - alpha) * (target_sd^2)
      term2 <- sum(1 / (1:i)) * tau
      mu_star <- sqrt(term1 + term2)
      
      data_list[[i]][,1] <- scaled_x + mu_star
      data_list[[i]][,2] <- scaled_y + mu_star
    }
  }
  
  # --- Flatten Data ---
  # Combine the normalized list back into a single N x 2 matrix for geometric calc
  flat_data <- do.call(rbind, data_list)
  
  # --- Validation & Core Calculation ---
  if (length(start) != 2) stop("Start must be a numeric vector of length 2")
  valid_methods <- c("min", "median", "equal", "mean", "neigh")
  if (!method %in% valid_methods) stop("Invalid method selected.")
  
  return(angle_alpha(flat_data, alpha, start, method))
}

#' Core angle calculation function
#' 
#' Internal function containing the geometric logic.
#' @keywords internal
angle_alpha <- function(data, alpha, start, method_select) {
  
  # Filter data > start
  data <- data[((data[,1] > start[1]) & (data[,2] > start[2])), , drop = FALSE]
  
  if (nrow(data) < 2) {
    stop("Not enough data points remaining after filtering (X > start[1] & Y > start[2]). Check your data or normalization parameters.")
  }
  
  outnum <- floor(dim(data)[1] * (1 - alpha))
  
  # --- Internal Helpers ---
  ab <- function(point, start){
    a <- (point[2]-start[2])/(point[1]-start[1])
    b <- (start[2]*point[1]-point[2]*start[1])/(point[1]-start[1])
    return(c(a,b))
  }
  
  abangle <- function(a,b,point,start){
    pointtrans <- point-start
    if(a>=0){
      angle <- atan(a)/pi*180
    }else if((a<0)&(pointtrans[1]<0)){
      angle <- atan(a)/pi*180+180
    }else{
      angle <- atan(a)/pi*180
    }
    return(angle)
  }
  
  # --- Pre-calculation Loop ---
  n <- nrow(data)
  Tan <- numeric(n)
  Angle <- numeric(n)
  Coef <- matrix(0, nrow=n, ncol=2)
  
  for(i in 1:n){
    pt <- data[i,]
    cf <- ab(pt, start)
    Tan[i] <- cf[1]
    Coef[i,] <- cf
    Angle[i] <- abangle(cf[1], cf[2], pt, start)
  }
  
  # --- Ordering & Pairing ---
  orderup <- order(-Angle)
  orderdown <- order(Angle)
  
  comple <- function(outnum, Angle, orderup, orderdown){
    updown <- matrix(0, nrow = outnum+1, ncol = 2)
    anglevs <- numeric(outnum+1)
    for(i in 1:(outnum+1)){
      up <- orderup[i]
      down <- orderdown[outnum+2-i]
      updown[i,] <- c(up, down)
      anglevs[i] <- abs(Angle[up]-Angle[down])
    }
    return(list(updown, anglevs))
  }
  
  # --- Logic Dispatch based on outnum ---
  angleminlist <- NULL; angleequallist <- NULL
  anglemeanlist <- NULL; anglemedianlist <- NULL; angleneighlist <- NULL
  
  if(outnum==0){
    whereup <- orderup[1]
    wheredown <- orderdown[1]
    anglemin <- abs(Angle[whereup]-Angle[wheredown])
    tanup <- Tan[whereup]; tandown <- Tan[wheredown]
    tan2 <- (tanup-tandown)/(1+tanup*tandown)
    
    row_res <- c(anglemin, tan2)
    res_mat <- rbind(row_res, data[whereup,], data[wheredown,], Coef[whereup,], Coef[wheredown,])
    
    # All methods return the same if outnum=0
    angleminlist <- res_mat; angleequallist <- res_mat
    anglemeanlist <- res_mat; anglemedianlist <- res_mat; angleneighlist <- res_mat
    
  } else {
    updownanglevs <- comple(outnum, Angle, orderup, orderdown)
    updown <- updownanglevs[[1]]
    anglevs <- updownanglevs[[2]]
    
    # 1. Min
    where <- order(anglevs)[1]
    whereup <- updown[where,1]; wheredown <- updown[where,2]
    anglemin <- anglevs[where]
    tanup <- Tan[whereup]; tandown <- Tan[wheredown]
    tan2 <- (tanup-tandown)/(1+tanup*tandown)
    angleminlist <- rbind(c(anglemin,tan2), data[whereup,], data[wheredown,], Coef[whereup,], Coef[wheredown,])
    
    # 2. Equal
    if(outnum==1){
      angleequallist <- angleminlist
    } else if(outnum%%2==0){
      where <- ceiling(outnum/2)
      whereup <- updown[where,1]; wheredown <- updown[where,2]
      angleequal <- anglevs[where]
      tanup <- Tan[whereup]; tandown <- Tan[wheredown]
      tan2 <- (tanup-tandown)/(1+tanup*tandown)
      angleequallist <- rbind(c(angleequal,tan2), data[whereup,], data[wheredown,], Coef[whereup,], Coef[wheredown,])
    } else {
      where <- floor(outnum/2)
      whereup <- updown[where,1]; whereupup <- updown[where+1,1]
      wheredown <- updown[where,2]; wheredowndown <- updown[where+1,2]
      
      pointup <- apply(data[c(whereup,whereupup),,drop=FALSE], 2, mean)
      coefup <- ab(pointup, start)
      pointdown <- apply(data[c(wheredown,wheredowndown),,drop=FALSE], 2, mean)
      coefdown <- ab(pointdown, start)
      
      angleequal <- abs(mean(Angle[c(whereup,whereupup)]) - mean(Angle[c(wheredown,wheredowndown)]))
      tan2 <- tan(angleequal*pi/180) # tan expects radians usually, but keeping original logic flow
      angleequallist <- rbind(c(angleequal,tan2), pointup, pointdown, coefup, coefdown)
    }
    
    # 3. Mean
    anglemean <- mean(anglevs)
    tan2 <- tan(anglemean*pi/180)
    pointup <- apply(data[updown[,1],,drop=FALSE], 2, mean)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[updown[,2],,drop=FALSE], 2, mean)
    coefdown <- ab(pointdown, start)
    anglemeanlist <- rbind(c(anglemean,tan2), pointup, pointdown, coefup, coefdown)
    
    # 4. Median
    anglemedian <- median(anglevs)
    tan2 <- tan(anglemedian*pi/180)
    pointup <- apply(data[updown[,1],,drop=FALSE], 2, median)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[updown[,2],,drop=FALSE], 2, median)
    coefdown <- ab(pointdown, start)
    anglemedianlist <- rbind(c(anglemedian,tan2), pointup, pointdown, coefup, coefdown)
    
    # 5. Neighbor
    where <- order(anglevs)[1]
    whereup <- updown[where,1]
    upnum <- which(orderup==whereup)
    whereupup <- orderup[ifelse(upnum==1, upnum, upnum-1)]
    
    wheredown <- updown[where,2]
    downnum <- which(orderdown==wheredown)
    wheredowndown <- orderdown[ifelse(downnum==1, downnum, downnum-1)]
    
    pointup <- apply(data[c(whereup,whereupup),,drop=FALSE], 2, mean)
    coefup <- ab(pointup, start)
    pointdown <- apply(data[c(wheredown,wheredowndown),,drop=FALSE], 2, mean)
    coefdown <- ab(pointdown, start)
    
    angleneigh <- abs(mean(Angle[c(whereup,whereupup)]) - mean(Angle[c(wheredown,wheredowndown)]))
    tan2 <- tan(angleneigh*pi/180)
    angleneighlist <- rbind(c(angleneigh,tan2), pointup, pointdown, coefup, coefdown)
  }
  
  # Method Mapping
  method_map <- c(
    "min" = "angleminlist",
    "median" = "anglemedianlist",
    "equal" = "angleequallist",
    "mean" = "anglemeanlist",
    "neigh" = "angleneighlist"
  )
  
  # Return variables from local environment
  result <- get(method_map[method_select])
  return(result)
}