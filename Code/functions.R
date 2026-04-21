# ==============================================================================
# File: functions.R
# Description: Defines Utility, Feasibility checks, and Objective Functions.
# ==============================================================================

check_parameters <- function(theta, omega) {
  stopifnot("Discount factor must be between 0 and 1" = theta$beta > 0 && theta$beta < 1)
  stopifnot("Risk aversion must be positive" = theta$gamma > 0)
  stopifnot("Minimum consumption must be positive" = theta$minCons > 0)
  stopifnot("Interest rate cannot be -1" = abs(1.0 + theta$r) > 1e-8)
}

min_and_max_assets <- function(theta, max_possible_income) {
  MinMaxA <- matrix(0, nrow = T_periods + 1, ncol = 2)
  MinMaxA[T_periods + 1, 1] <- 0.0 # Cannot die with debt
  
  for (i in T_periods:1) {
    MinMaxA[i, 1] <- MinMaxA[i + 1, 1] / (1.0 + theta$r) + theta$minCons
  }
  
  MinMaxA[1, 2] <- theta$A0 + max_possible_income
  for (t in 2:(T_periods + 1)) {
    MinMaxA[t, 2] <- (MinMaxA[t - 1, 2] + max_possible_income - theta$minCons) * (1.0 + theta$r)
  }
  return(MinMaxA)
}

u <- function(c, theta) {
  if (c < theta$minCons) c <- theta$minCons
  if (abs(theta$gamma - 1.0) < 1e-8) {
    return(log(c))
  } else {
    return((c^(1.0 - theta$gamma)) / (1.0 - theta$gamma))
  }
}

#' Objective Function for Workers (Faces Markov Income Shocks)
obj_func_worker <- function(A1, A0, IncomeNow, Kchoice, VNext, TM, y, theta) {
  c <- A0 + IncomeNow - A1 / (1.0 + theta$r)
  V_tomorrow <- 0.0
  for (i in 1:NY) {
    V_tomorrow <- V_tomorrow + TM[y, i] * interp1D(Kchoice, VNext[, i], A1)
  }
  return( -(u(c, theta) + theta$beta * V_tomorrow) )
}

#' Objective Function for Students (Deterministic income, Continuation Value is 1D)
obj_func_student <- function(A1, A0, IncomeNow, Kchoice, VNext1D, theta) {
  c <- A0 + IncomeNow - A1 / (1.0 + theta$r)
  V_tomorrow <- interp1D(Kchoice, VNext1D, A1)
  return( -(u(c, theta) + theta$beta * V_tomorrow) )
}