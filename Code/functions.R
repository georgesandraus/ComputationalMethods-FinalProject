# ==============================================================================
# File: functions.R
# Description: Defines Utility, Feasibility checks, and the Objective Function.
# ==============================================================================

#' Checks if parameters are economically and mathematically feasible
check_parameters <- function(theta, omega) {
  stopifnot("Discount factor must be between 0 and 1" = theta$beta > 0 && theta$beta < 1)
  stopifnot("Risk aversion must be positive" = theta$gamma > 0)
  stopifnot("Minimum consumption must be positive" = theta$minCons > 0)
  stopifnot("Interest rate cannot be -1" = abs(1.0 + theta$r) > 1e-8)
}

#' Defines the minimum and maximum assets for the grid
#' @param theta Parameter list
#' @param max_possible_income The highest possible income in the grid
#' @return A matrix with bounds for each period
min_and_max_assets <- function(theta, max_possible_income) {
  MinMaxA <- matrix(0, nrow = T_periods + 1, ncol = 2)
  
  # 1. Borrowing Constraint (Minimum Assets)
  MinMaxA[T_periods + 1, 1] <- 0.0 # Cannot die with debt
  
  for (i in T_periods:1) {
    # Strictest borrowing constraint: Must afford minCons even if income is zero
    MinMaxA[i, 1] <- MinMaxA[i + 1, 1] / (1.0 + theta$r) + theta$minCons
  }
  
  # 2. Maximum Assets (Assuming saving max income and starving)
  MinMaxA[1, 2] <- theta$A0 + max_possible_income
  
  for (t in 2:(T_periods + 1)) {
    MinMaxA[t, 2] <- (MinMaxA[t - 1, 2] + max_possible_income - theta$minCons) * (1.0 + theta$r)
  }
  
  return(MinMaxA)
}

#' Per-period Utility Function (CRRA)
u <- function(c, theta) {
  # Defensive: Penalize heavily if consumption drops below minimum
  if (c < theta$minCons) {
    c <- theta$minCons
  }
  
  if (abs(theta$gamma - 1.0) < 1e-8) {
    return(log(c))
  } else {
    return((c^(1.0 - theta$gamma)) / (1.0 - theta$gamma))
  }
}

#' Objective Function for the Optimizer (Inside the Bellman Equation)
#' Returns the NEGATIVE expected value (because R's optimize() finds minimums)
#' @param A1 Chosen savings for tomorrow
#' @param A0 Current assets today
#' @param IncomeNow Current realized income today
#' @param Kchoice Grid of possible assets tomorrow
#' @param VNext Value function matrix for tomorrow (dim: NK x NY)
#' @param TM Transition Matrix for the specific income process (dim: NY x NY)
#' @param y Current income state index (1 to NY)
obj_func <- function(A1, A0, IncomeNow, Kchoice, VNext, TM, y, theta) {
  
  # Calculate consumption implied by the budget constraint
  c <- A0 + IncomeNow - A1 / (1.0 + theta$r)
  
  # Calculate the Expected Value of tomorrow
  V_tomorrow <- 0.0
  for (i in 1:NY) {
    # Interpolate VNext for the chosen saving A1 in state i
    V_interp <- interp1D(Kchoice, VNext[, i], A1)
    # Multiply by the transition probability from state y to state i
    V_tomorrow <- V_tomorrow + TM[y, i] * V_interp
  }
  
  # Return negative because optimize() minimizes by default
  return( -(u(c, theta) + theta$beta * V_tomorrow) )
}