# ==============================================================================
# File: parameters.R
# Description: Defines the economic environment, parameters, and grids.
# ==============================================================================

# Environment 
T_periods <- 50    # Number of time periods (e.g., age 15 to 65)
T_school  <- 3     # Periods where studying is a choice (t = 1, 2, 3)
NK        <- 500    # Number of points in the asset grid
NY        <- 3     # Number of states in the income Markov process
NSIM      <- 10000  # Number of simulated individuals

# Preferences and Financial Parameters (theta)
theta <- list(
  beta    = 0.95,  # Discount factor (could be lowered later to test myopia)
  gamma   = 1.5,   # Coefficient of relative risk aversion
  r       = 0.05,  # Interest rate
  A0      = 0.5,   # Initial assets (low wealth scenario)
  minCons = 1e-4   # Minimum consumption (to avoid log(0) or division by zero)
)

# Income, Policy, and Markov Parameters (omega)
omega <- list(
  # Base wages
  w_L = 1.0,       # Base wage for Low-Skilled (in minimum wages)
  w_H = 2.5,       # Base wage for High-Skilled
  
  # Education parameters
  cost_edu = 0.3,  # Opportunity/direct cost of studying per period
  benefit  = 0.2,  # "Pé-de-Meia" conditional cash transfer per period, with
                   # an approx. value of .2 minimum wages
  
  # Markov Process for Low-Skilled (High Variance / High Risk)
  Y_L_grid = c(0.6, 1.0, 1.4),  # Bad, Normal, Good states
  TM_L = matrix(c(
    0.4, 0.4, 0.2,  # From Bad
    0.2, 0.6, 0.2,  # From Normal
    0.2, 0.4, 0.4   # From Good
  ), nrow = 3, byrow = TRUE),
  
  # Markov Process for High-Skilled (Low Variance / High Stability)
  Y_H_grid = c(2.2, 2.5, 2.8),
  TM_H = matrix(c(
    0.6, 0.3, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.3, 0.6
  ), nrow = 3, byrow = TRUE),
  
  Yint_0 = 2 # Everyone starts at the 'Normal' state (index 2)
)

# Computational Settings
grid_spacing  <- "logsteps" # "equalsteps" or "logsteps"
interp_method <- "linear"   # We'll use linear via 'approx' in R for robust monotonic interpolation