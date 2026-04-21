# ==============================================================================
# File: solveModel.R
# Description: Solves the model using Backward Induction. Handles the sequential
#              dropout option for students.
# ==============================================================================

# Solves the standard working life-cycle (Low-Skilled or High-Skilled)
solve_worker_path <- function(theta, omega, Kgrid, skill) {
  if (skill == "L") {
    Income_Grid <- omega$Y_L_grid; TM <- omega$TM_L; wage_base <- omega$w_L
    start_t <- 1
  } else {
    Income_Grid <- omega$Y_H_grid; TM <- omega$TM_H; wage_base <- omega$w_H
    start_t <- T_school + 1 # High-Skilled only starts at t=4
  }
  
  V <- array(0, dim = c(T_periods + 1, NK, NY))
  gK <- array(0, dim = c(T_periods, NK, NY))
  
  for (age in T_periods:start_t) {
    for (k in 1:NK) {
      for (y in 1:NY) {
        IncomeNow <- wage_base * Income_Grid[y]
        AssetsNow <- Kgrid[age, k]
        
        minK <- Kgrid[age + 1, 1]
        maxK <- (AssetsNow + IncomeNow - theta$minCons) * (1.0 + theta$r)
        
        if (maxK <= minK) {
          gK[age, k, y] <- minK; V[age, k, y] <- -1e6
          next
        }
        
        res <- optimize(f = obj_func_worker, interval = c(minK, maxK),
                        A0 = AssetsNow, IncomeNow = IncomeNow, Kchoice = Kgrid[age + 1, ],
                        VNext = V[age + 1, , ], TM = TM, y = y, theta = theta)
        gK[age, k, y] <- res$minimum
        V[age, k, y]  <- -res$objective
      }
    }
  }
  return(list(V = V, gK = gK))
}

# Solves the 3-year schooling phase, incorporating the sequential option to drop out
solve_student_path <- function(theta, omega, Kgrid, V_W, V_H) {
  V_S   <- matrix(0, nrow = T_school, ncol = NK)
  gK_S  <- matrix(0, nrow = T_school, ncol = NK)
  V_Opt <- matrix(0, nrow = T_school, ncol = NK) # The Option Value (Study or Work)
  
  costs_t123 <- omega$benefit_path - omega$cost_edu
  
  for (age in T_school:1) {
    IncomeNow <- costs_t123[age]
    
    # What happens tomorrow?
    if (age == T_school) {
      # If they finish t=3, they enter High-Skilled market at t=4 in 'Normal' state
      VNext1D <- V_H[age + 1, , omega$Yint_0]
    } else {
      # If they finish t=1 or t=2, they face the option to stay or drop out at t+1
      VNext1D <- V_Opt[age + 1, ]
    }
    
    for (k in 1:NK) {
      AssetsNow <- Kgrid[age, k]
      minK <- Kgrid[age + 1, 1]
      maxK <- (AssetsNow + IncomeNow - theta$minCons) * (1.0 + theta$r)
      
      if (maxK <= minK) {
        gK_S[age, k] <- minK; V_S[age, k] <- -1e6
      } else {
        res <- optimize(f = obj_func_student, interval = c(minK, maxK),
                        A0 = AssetsNow, IncomeNow = IncomeNow, 
                        Kchoice = Kgrid[age + 1, ], VNext1D = VNext1D, theta = theta)
        gK_S[age, k] <- res$minimum
        V_S[age, k]  <- -res$objective
      }
      
      # The Real Option: Compare value of studying this year vs dropping out immediately
      V_W_dropout   <- V_W[age, k, omega$Yint_0]
      V_Opt[age, k] <- max(V_S[age, k], V_W_dropout)
    }
  }
  return(list(V_S = V_S, gK_S = gK_S, V_Opt = V_Opt))
}