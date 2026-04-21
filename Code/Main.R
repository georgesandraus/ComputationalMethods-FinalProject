# ==============================================================================
# Computational Methods in Economics - Final Project
# Author: Georges Mikhael Andraus
# Description: Evaluating the "Pé-de-Meia" Program via a Life-Cycle Model 
#              with Sequential Discrete Choices and Stochastic Income.
# ==============================================================================

# ==============================================================================
# 1. Setup and Environment Management
# ==============================================================================
# Package management via 'renv'
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
# renv::activate() 

library(tidyverse)
library(patchwork)
library(microbenchmark)
library(stargazer)

# Set seed for reproducibility
set.seed(12345)

# Ensure output directories exist
dir.create("Output/Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Output/Tables", showWarnings = FALSE, recursive = TRUE)

# Load external modules
source("Code/parameters.R")
source("Code/procedures.R")
source("Code/functions.R")
source("Code/solveModel.R")

# ==============================================================================
# 2. Initialization and Sanity Checks
# ==============================================================================
cat("\n--- Phase 1: Initializing Model Constraints ---\n")

check_parameters(theta, omega)

max_possible_income <- max(
  omega$w_L * max(omega$Y_L_grid), 
  omega$w_H * max(omega$Y_H_grid)
)

MinMaxA <- min_and_max_assets(theta, max_possible_income)
Kgrid   <- get_K_Grid(theta, MinMaxA)

cat("Asset Grid dimensions:", nrow(Kgrid), "periods x", ncol(Kgrid), "nodes.\n")

# ==============================================================================
# 3. Model Solution (Backward Induction for Sequential Paths)
# ==============================================================================
cat("\n--- Phase 2: Solving Dynamic Programming Paths ---\n")

# Path 1: Worker (Low-Skilled) from t=1 to T
cat("Solving Worker Life-Cycle (Low-Skilled)...\n")
sol_W <- solve_worker_path(theta, omega, Kgrid, skill = "L")

# Path 2: Worker (High-Skilled) from t=4 to T
cat("Solving Worker Life-Cycle (High-Skilled)...\n")
sol_H <- solve_worker_path(theta, omega, Kgrid, skill = "H")

# Path 3: Student (t=1, 2, 3) with Sequential Dropout Option
cat("Solving Student Path (Sequential Dropout Option)...\n")
sol_S <- solve_student_path(theta, omega, Kgrid, sol_W$V, sol_H$V)

# ==============================================================================
# 4. The Initial Discrete Choice at t=1
# ==============================================================================
cat("\n--- Phase 3: Evaluating Initial Discrete Choice ---\n")

V1_W <- sol_W$V[1, , omega$Yint_0]
V1_S <- sol_S$V_S[1, ]

df_choice <- data.frame(
  Initial_Wealth = Kgrid[1, ],
  V_Work = V1_W,
  V_Study = V1_S
) %>%
  filter(V_Work > -1e5 & V_Study > -1e5) %>%
  mutate(Optimal_Choice = ifelse(V_Study >= V_Work, "Study", "Work"))

threshold_idx <- which(df_choice$V_Study >= df_choice$V_Work)[1]
A_star <- ifelse(is.na(threshold_idx), "Never", 
                 round(df_choice$Initial_Wealth[threshold_idx], 3))

cat("Minimum Initial Wealth required to study (Threshold A*):", A_star, "\n")

p_choice <- ggplot(df_choice, aes(x = Initial_Wealth)) +
  geom_line(aes(y = V_Work, color = "Work (Dropout)"), linewidth = 1.2) +
  geom_line(aes(y = V_Study, color = "Study"), linewidth = 1.2) +
  geom_vline(xintercept = as.numeric(A_star), linetype = "dashed", 
             color = "black") +
  annotate("text", x = as.numeric(A_star) + 0.5, y = -15, 
           label = paste("Threshold A* =", A_star), fontface = "bold") +
  coord_cartesian(ylim = c(-60, 0), xlim = c(0, 3.5)) + 
  labs(
    title = "",
    subtitle = "",
    x = "Initial Wealth (Assets)",
    y = "Value Function V(a, t=1)",
    color = "Choice"
  ) +
  scale_color_manual(values = c("Study" = "#E64B35", 
                                "Work (Dropout)" = "#4DBBD5")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Output/Figures/discrete_choice_t1.pdf", p_choice, width = 8, height = 5, dpi = 300)

# ==============================================================================
# 5. Computational Benchmarking (N_Grid vs Execution Time)
# ==============================================================================
cat("\n--- Phase 4: Computational Benchmarking ---\n")
cat("Running Microbenchmark (NK=100 vs NK=500). Please wait, this takes time...\n")

benchmark_model <- function(test_NK) {
  assign("NK", test_NK, envir = .GlobalEnv)
  test_Kgrid <- get_K_Grid(theta, MinMaxA)
  invisible(solve_worker_path(theta, omega, test_Kgrid, skill = "L"))
}

mb_results <- microbenchmark(
  "NK = 100" = benchmark_model(100),
  "NK = 500" = benchmark_model(500),
  times = 1 
)

assign("NK", 500, envir = .GlobalEnv) # Restore high precision grid

df_bench <- summary(mb_results, unit = "s") %>%
  mutate(
    Grid_Size = expr,
    Execution_Time_Seconds = round(mean, 2)
  ) %>%
  select(Grid_Size, Execution_Time_Seconds)

print(df_bench)

stargazer(df_bench, type = "latex", summary = FALSE, rownames = FALSE,
          title = "Computational Performance by Asset Grid Size",
          out = "Output/Tables/benchmark_results.tex", float = FALSE)


# ==============================================================================
# 6. Mechanism: Monte Carlo Simulation of Paths (Sequential Dropout)
# ==============================================================================
cat("\n--- Phase 5: Simulating Life-Cycle Paths with Dropout Mechanism ---\n")

simulate_agent <- function(initial_wealth) {
  a <- numeric(T_periods + 1); c <- numeric(T_periods); status <- rep("W", 
                                                                      T_periods)
  a[1] <- initial_wealth
  currently_studying <- TRUE
  
  for (t in 1:T_periods) {
    # Check Sequential Dropout Option
    if (currently_studying && t <= T_school) {
      v_s <- interp1D(Kgrid[t, ], sol_S$V_S[t, ], a[t])
      v_w <- interp1D(Kgrid[t, ], sol_W$V[t, , omega$Yint_0], a[t])
      if (v_w > v_s) currently_studying <- FALSE
    }
    
    if (currently_studying && t <= T_school) {
      status[t] <- "S"
      inc <- omega$benefit_path[t] - omega$cost_edu
      a1 <- interp1D(Kgrid[t+1, ], sol_S$gK_S[t, ], a[t])
    } else if (currently_studying && t > T_school) {
      status[t] <- "H"
      inc <- omega$w_H * omega$Y_H_grid[omega$Yint_0]
      a1 <- interp1D(Kgrid[t+1, ], sol_H$gK[t, , omega$Yint_0], a[t])
    } else {
      status[t] <- "L"
      inc <- omega$w_L * omega$Y_L_grid[omega$Yint_0]
      a1 <- interp1D(Kgrid[t+1, ], sol_W$gK[t, , omega$Yint_0], a[t])
    }
    
    a[t+1] <- a1
    c[t] <- a[t] + inc - a1 / (1.0 + theta$r)
  }
  return(data.frame(Age = 1:T_periods, Consumption = c, 
                    Assets = a[1:T_periods], 
                    Status = status, A0 = initial_wealth))
}

sim_rich <- simulate_agent(1.0)
sim_poor <- simulate_agent(0.3)

p_mechanisms <- bind_rows(sim_rich, sim_poor) %>%
  mutate(Profile = paste("A0 =", A0, "| Ends as:", Status[T_periods])) %>%
  pivot_longer(cols = c(Consumption, Assets), 
               names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Age, y = Value, color = Profile, linetype = Profile)) +
  geom_line(linewidth = 1) + facet_wrap(~Variable, scales = "free_y", ncol=1) +
  labs(title = "", 
       x = "Age (Periods)") + theme_minimal()

ggsave("Output/Figures/mechanism_paths.pdf", p_mechanisms, 
       width = 8, height = 7, dpi = 300)


# ==============================================================================
# 7. Policy Counterfactuals: The Timing of Benefits
# ==============================================================================
cat("\n--- Phase 6: Testing Different Program Designs (Timing) ---\n")

policy_designs <- list(
  "No Policy"             = c(0.0, 0.0, 0.0),
  "Even (Baseline)"       = c(0.2, 0.2, 0.2),
  "End-Loaded (Bonus)"    = c(0.0, 0.0, 0.6), # Big bonus at graduation
  "Front-Loaded (Rescue)" = c(0.4, 0.1, 0.1)  # Heavy support in the 1st year
)

df_policy_V <- data.frame()
thresholds <- list()

for (p_name in names(policy_designs)) {
  temp_omega <- omega
  temp_omega$benefit_path <- policy_designs[[p_name]]
  
  # Re-solve ONLY the Study path (Work paths sol_W and sol_H are reused!)
  t_sol_S <- solve_student_path(theta, temp_omega, Kgrid, sol_W$V, sol_H$V)
  V_S_temp <- t_sol_S$V_S[1, ]
  
  diff <- V_S_temp - V1_W
  idx <- which(diff >= 0)[1]
  A_star <- ifelse(is.na(idx), NA, round(Kgrid[1, idx], 3))
  thresholds[[p_name]] <- A_star
  
  df_policy_V <- bind_rows(df_policy_V, data.frame(
    Initial_Wealth = Kgrid[1, ], V_Study = V_S_temp, Policy = p_name
  ))
}

df_policy_results <- data.frame(
  Program_Design = names(thresholds),
  Threshold_A_star = unlist(thresholds)
)
print(df_policy_results)

stargazer(df_policy_results, type = "latex", summary = FALSE, rownames = FALSE,
          title = "Policy Design Comparison: Minimum Wealth Required to Study",
          out = "Output/Tables/policy_designs.tex", float = FALSE)

# Plotting the Value Functions of the Policies
df_plot_policy <- df_policy_V %>% filter(Initial_Wealth <= 2.0 & V_Study > -100)
df_plot_W <- df_choice %>% filter(Initial_Wealth <= 2.0 & V_Work > -100)

p_designs <- ggplot() +
  geom_line(data = df_plot_W, aes(x = Initial_Wealth, y = V_Work, 
                                  linetype = "Work (Outside Option)"), 
            color = "black", linewidth = 1) +
  geom_line(data = df_plot_policy, aes(x = Initial_Wealth, y = V_Study, 
                                       color = Policy), linewidth = 1.2) +
  coord_cartesian(ylim = c(-30, -5)) +
  labs(
    title = "Impact of Benefit Timing on Lifetime Utility",
    subtitle = "Front-loading helps the poorest cross the threshold; End-loading fails due to credit constraints.",
    x = "Initial Wealth (Assets)", y = "Expected Lifetime Utility at t=1", color = "Program Design", linetype = ""
  ) +
  theme_minimal() + theme(legend.position = "right")

ggsave("Output/Figures/policy_designs_comparison.pdf", p_designs, width = 9, height = 5, dpi = 300)

cat("\nProject Code Execution Completed Successfully. Outputs saved to /Output.\n")