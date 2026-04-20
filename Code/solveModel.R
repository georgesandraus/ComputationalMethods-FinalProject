# ==============================================================================
# File: solveModel.R
# Description: Solves the dynamic programming problem using Backward Induction.
#              Computes V(a) for both the Work (L) and Study (H) paths.
# ==============================================================================

#' Solves the finite-horizon consumption-savings model
#' @param theta Economic parameters
#' @param omega Income and Markov parameters
#' @param Kgrid The asset grid
#' @param path String: "W" (Work/Low-Skilled) or "S" (Study/High-Skilled)
#' @return A list containing the Value function (V) and Policy function (gK)
solve_model_path <- function(theta, omega, Kgrid, path) {
  
  # 1. Setup Path-Specific Variables
  if (path == "W") {
    Income_Grid <- omega$Y_L_grid
    TM          <- omega$TM_L
    wage_base   <- omega$w_L
    # No costs or benefits for the working path
    costs_t123  <- rep(0.0, 3) 
    
  } else if (path == "S") {
    Income_Grid <- omega$Y_H_grid
    TM          <- omega$TM_H
    wage_base   <- omega$w_H
    # Net income in t=1,2,3 for students: (Benefit - Cost). 
    # Negative value means they must borrow or rely on initial wealth.
    costs_t123  <- rep(omega$benefit - omega$cost_edu, 3) 
  } else {
    stop("Path must be 'W' or 'S'")
  }
  
  # Initialize Value and Policy Matrices
  # Dimensions: Time x Asset Grid x Income State
  V  <- array(0, dim = c(T_periods + 1, NK, NY))
  gK <- array(0, dim = c(T_periods, NK, NY))
  
  # 2. Terminal Condition (Period T+1)
  # No utility from leaving assets behind (unless bequest motive is added)
  V[T_periods + 1, , ] <- 0.0 
  
  # 3. Backward Induction Loop
  for (age in T_periods:1) {
    
    # In t=1,2,3, the Student path has fixed (negative) income and NO Markov shocks yet.
    # The Worker path has standard Markov shocks from day 1.
    is_schooling_period <- (path == "S" && age <= T_school)
    
    for (k in 1:NK) {
      for (y in 1:NY) {
        
        # Determine current income
        if (is_schooling_period) {
          # While in school, income is deterministically the net subsidy
          IncomeNow <- costs_t123[age]
          # Transition matrix is technically irrelevant during school, 
          # but we need to compute expected value of entering the labor market in t=4
          # Assuming they enter the market at the 'Normal' state (y=2)
          effective_TM <- matrix(0, nrow=3, ncol=3)
          effective_TM[, omega$Yint_0] <- 1.0 
        } else {
          # After school (or always for workers), income is subject to Markov shocks
          IncomeNow <- wage_base * Income_Grid[y]
          effective_TM <- TM
        }
        
        # Current Assets and Choice Grid
        AssetsNow <- Kgrid[age, k]
        Kchoice   <- Kgrid[age + 1, ]
        VNext     <- V[age + 1, , ]
        
        # Borrowing Constraint limits
        minK <- Kchoice[1]
        
        # Max feasible savings: Save all current resources minus minimum consumption
        maxK <- (AssetsNow + IncomeNow - theta$minCons) * (1.0 + theta$r)
        
        # Defensive Check: If resources are so low that they can't even afford minCons
        if (maxK <= minK) {
          gK[age, k, y] <- minK
          # Penalty for starvation (large negative utility)
          V[age, k, y]  <- -1e6
          next
        }
        
        # 4. Optimization (Brent's Method via optimize)
        # optimize() searches for a minimum. obj_func returns negative expected utility.
        res <- optimize(
          f = obj_func,
          interval = c(minK, maxK),
          A0 = AssetsNow,
          IncomeNow = IncomeNow,
          Kchoice = Kchoice,
          VNext = VNext,
          TM = effective_TM,
          y = y,
          theta = theta
        )
        
        # 5. Store optimal policy and resulting value
        gK[age, k, y] <- res$minimum
        V[age, k, y]  <- -res$objective # Convert back to positive utility
        
      }
    }
  }
  
  return(list(V = V, gK = gK))
}