# Evaluating the "Pé-de-Meia" Program: A Life-Cycle Model with Sequential Discrete Choices

**Author:** Georges Mikhael Andraus, Sao Paulo School of Economics - FGV

---

## 1. Overview

This repository contains the replication package for the final project of the **Computational Methods in Economics** course. 

The project evaluates the impact of the Brazilian "Pé-de-Meia" program on educational attainment through a structural life-cycle model. The framework features a sequential discrete choice (stay in school vs. drop out) over the first three periods of life, subject to stochastic income shocks and strict borrowing constraints. The goal is to identify how the timing and magnitude of financial incentives alter the optimal human capital investment decisions of credit-constrained agents.

---

## 2. Repository Structure

```text
.
├── Code/
│   ├── parameters.R            # Defines model environment, economic params, and grids
│   ├── functions.R             # Core utilities: CRRA utility and path-specific objective functions
│   ├── procedures.R            # Routines for asset grid construction and robust interpolation
│   ├── solveModel.R            # Main solver implementing Backward Induction and Real Options logic
│   └── Main.R                  # Master script: runs benchmarks, simulations, and counterfactuals
├── Output/
│   ├── Figures/                # Generated plots (Discrete choice, Mechanisms, Policy designs)
│   └── Tables/                 # Generated LaTeX tables (Benchmarks, Policy results)
├── renv/                       # R environment isolation files
├── renv.lock                   # Exact package versions for reproducibility
├── ComputationalMethods-FinalProject.Rproj 
└── README.md
```

### Code Descriptions

*   **`parameters.R`**: Centralizes all economic primitives (wages, discount factors, costs) and computational settings (grid size $N_K=500$).
*   **`functions.R`**: Contains the CRRA utility function and the specific objective functions for workers and students used by the optimizer.
*   **`procedures.R`**: Handles the construction of log-spaced grids to capture curvature near the borrowing constraint and robust linear interpolation.
*   **`solveModel.R`**: The core computational engine. It solves the life-cycle problem via Backward Induction. For the schooling phase ($t \le 3$), it incorporates the sequential dropout option by evaluating the value of the "Real Option" at each state.
*   **`Main.R`**: The driver script. It executes the full model, runs the computational benchmark (performance testing), performs policy counterfactuals, and exports all visual and tabular outputs.

---

## 3. Computational Environment

The analysis was conducted using the following environment:
*   **Operating System:** macOS Tahoe 26.4.1
*   **R Version:** 4.5.1
*   **RStudio Version:** 2026.01.1+403

---

## 4. Replication Instructions

This project uses `renv` to ensure package version consistency. To replicate the results:

1.  Open the `ComputationalMethods-FinalProject.Rproj` file in RStudio.
2.  Install the `renv` package if not already present: `install.packages("renv")`.
3.  Run **`renv::restore()`** in the R console. This will automatically install the exact versions of `tidyverse`, `stargazer`, `microbenchmark`, and `patchwork` used in this project.
4.  Run the **`Code/Main.R`** script. 
    *   *Note:* Given the high-precision grid ($N_K=500$), the full execution including benchmarking and counterfactuals may take a different time depending on your CPU.

**AI Usage Declaration:**  
Portions of the R code and the LaTeX structuring of this project were developed with the assistance of an LLM (Gemini 3.1 Pro Preview). The AI was primarily used to enforce good computational practices (environment management via `renv`, functional programming, and optimizing grid evaluations). All economic modeling choices, parameter calibration, and analytical interpretations remain my own.
