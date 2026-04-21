# ==============================================================================
# File: procedures.R
# Description: Helper routines for Grid construction and Interpolation.
# ==============================================================================

# Generates the Asset Grid (KGrid)
# Argument: theta - List of economic parameters
# Argument: MinMaxA - Matrix (T_periods+1 x 2) with minimum and maximum feasible assets
# Function: A matrix of dimensions (T_periods+1 x NK)
get_K_Grid <- function(theta, MinMaxA) {
  
  Kgrid <- matrix(0, nrow = T_periods + 1, ncol = NK)
  Klog  <- matrix(0, nrow = T_periods + 1, ncol = NK)
  
  Kgrid[, 1]  <- MinMaxA[, 1]
  Kgrid[, NK] <- MinMaxA[, 2]
  
  Span <- Kgrid[, NK] - Kgrid[, 1]
  
  for (i in 1:(T_periods + 1)) {
    Klog[i, 1]  <- 0.0
    Klog[i, NK] <- Span[i]
  }
  
  # Build grid spacing
  for (i in 1:(T_periods + 1)) {
    for (j in 2:NK) {
      if (grid_spacing == "equalsteps") {
        Klog[i, j] <- Klog[i, j-1] + Span[i] / (NK - 1)
      } else if (grid_spacing == "logsteps") {
        Klog[i, j] <- Klog[i, j-1] + log(Span[i] + 1.0) / (NK - 1)
      } else {
        stop("grid_spacing must be 'equalsteps' or 'logsteps'")
      }
    }
  }
  
  # Transform back to levels
  for (i in 1:(T_periods + 1)) {
    for (j in 1:NK) {
      if (grid_spacing == "equalsteps") {
        Kgrid[i, j] <- Klog[i, j] + Kgrid[i, 1]
      } else if (grid_spacing == "logsteps") {
        Kgrid[i, j] <- exp(Klog[i, j]) - 1.0 + Kgrid[i, 1]
      }
    }
  }
  
  return(Kgrid)
}

# 1D Interpolation
# Argument: Xgrid - Vector of x values (asset grid)
# Argument: Ygrid - Vector of y values (value function)
# Argument: x - Target point to interpolate
# Function: Interpolated value
interp1D <- function(Xgrid, Ygrid, x) {
  # Defensive check: ensure no NAs in interpolation data
  if(any(is.na(Ygrid))) stop("NAs found in Ygrid during interpolation")
  
  if (interp_method == "linear") {
    # rule = 2 ensures that if the optimizer lies slightly out of bounds, 
    # it returns the nearest extreme value without throwing an NA
    res <- approx(x = Xgrid, y = Ygrid, xout = x, method = "linear", rule = 2)$y
  } else {
    stop("Interpolation method not supported.")
  }
  return(res)
}