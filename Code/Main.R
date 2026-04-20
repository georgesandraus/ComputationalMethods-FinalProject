# ==============================================================================
# Computational Methods in Economics - Final Project
# Author: Georges Mikhael Andraus
# Description: Evaluating the "Pé-de-Meia" Program via a Life-Cycle Model 
#              with Initial Discrete Choice and Stochastic Income.
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

# Check parameter feasibility (e.g. positive risk aversion, valid beta)
check_parameters(theta, omega)

# Determine the absolute maximum possible income across all paths to set the Asset Grid upper bound
max_possible_income <- max(
  omega$w_L * max(omega$Y_L_grid), 
  omega$w_H * max(omega$Y_H_grid)
)

# Get Min and Max Asset bounds
MinMaxA <- min_and_max_assets(theta, max_possible_income)

# Build the Log-spaced Asset Grid
Kgrid <- get_K_Grid(theta, MinMaxA)

cat("Asset Grid dimensions:", nrow(Kgrid), "periods x", ncol(Kgrid), "nodes.\n")

# ==============================================================================
# 3. Model Solution (Backward Induction for both Paths)
# ==============================================================================
cat("\n--- Phase 2: Solving Dynamic Programming Paths ---\n")

# Solve the entire life-cycle assuming the agent enters the labor market immediately (Low-Skilled)
cat("Solving Work Path (W)...\n")
sol_W <- solve_model_path(theta, omega, Kgrid, path = "W")

# Solve the entire life-cycle assuming the agent studies for t=1,2,3 (High-Skilled)
cat("Solving Study Path (S)...\n")
sol_S <- solve_model_path(theta, omega, Kgrid, path = "S")

# ==============================================================================
# 4. The Discrete Choice at t=1
# ==============================================================================
cat("\n--- Phase 3: Evaluating Initial Discrete Choice ---\n")

V1_W <- sol_W$V[1, , omega$Yint_0]
V1_S <- sol_S$V[1, , omega$Yint_0]

df_choice <- data.frame(
  Initial_Wealth = Kgrid[1, ],
  V_Work = V1_W,
  V_Study = V1_S
) %>%
  filter(V_Work > -1e5 & V_Study > -1e5) %>%
  mutate(Optimal_Choice = ifelse(V_Study >= V_Work, "Study", "Work"))

threshold_idx <- which(df_choice$V_Study >= df_choice$V_Work)[1]
A_star <- ifelse(is.na(threshold_idx), "Never", round(df_choice$Initial_Wealth[threshold_idx], 3))

cat("Minimum Initial Wealth required to study (Threshold A*):", A_star, "\n")

# Plotting the Value Functions at t=1 (WITH ZOOM)
p_choice <- ggplot(df_choice, aes(x = Initial_Wealth)) +
  geom_line(aes(y = V_Work, color = "Work (Low-Skilled)"), linewidth = 1.2) +
  geom_line(aes(y = V_Study, color = "Study (High-Skilled)"), linewidth = 1.2) +
  geom_vline(xintercept = as.numeric(A_star), linetype = "dashed", color = "black") +
  # Place the annotation dynamically near the threshold
  annotate("text", x = as.numeric(A_star) + 0.5, y = -15, 
           label = paste("Threshold A* =", A_star), fontface = "bold") +
  # ZOOM IN on the relevant area to avoid the starvation drop distortion
  coord_cartesian(ylim = c(-60, 0), xlim = c(0, 3.5)) + 
  labs(
    title = "Expected Lifetime Utility at t=1",
    subtitle = "Zoomed in. The vertical line is the minimum wealth needed to afford education.",
    x = "Initial Wealth (Assets)",
    y = "Value Function V(a, t=1)",
    color = "Choice"
  ) +
  scale_color_manual(values = c("Study (High-Skilled)" = "#E64B35", "Work (Low-Skilled)" = "#4DBBD5")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Output/Figures/discrete_choice_t1.png", p_choice, width = 8, height = 5, dpi = 300)

# ==============================================================================
# 5. Computational Benchmarking (N_Grid vs Execution Time)
# ==============================================================================
cat("\n--- Phase 4: Computational Benchmarking ---\n")
cat("Running Microbenchmark for different Grid Sizes (NK)...\n")

benchmark_model <- function(test_NK) {
  assign("NK", test_NK, envir = .GlobalEnv)
  test_Kgrid <- get_K_Grid(theta, MinMaxA)
  invisible(solve_model_path(theta, omega, test_Kgrid, path = "W"))
}

mb_results <- microbenchmark(
  "NK = 30" = benchmark_model(30),
  "NK = 60" = benchmark_model(60),
  times = 3 
)

assign("NK", 60, envir = .GlobalEnv)
print(mb_results)

# Fix: Force summary to use seconds so we don't crush the numbers
df_bench <- summary(mb_results, unit = "s") %>%
  mutate(
    Grid_Size = expr,
    Mean_Time_Seconds = round(mean, 2),
    Min_Time_Seconds  = round(min, 2),
    Max_Time_Seconds  = round(max, 2)
  ) %>%
  select(Grid_Size, Mean_Time_Seconds, Min_Time_Seconds, Max_Time_Seconds)

stargazer(df_bench, type = "latex", summary = FALSE, rownames = FALSE,
          title = "Computational Performance: Execution Time by Asset Grid Size",
          out = "Output/Tables/benchmark_results.tex",
          float = FALSE)

cat("\nProject Code Execution Completed Successfully. Outputs saved to /Output.\n")