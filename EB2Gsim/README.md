# Efford & Boulanger Two-Group Simulation

This folder extends the Efford & Boulanger (2019) simulation to a two-group setting with fixed grid spacing at `2 * sigma1`. The workflow keeps the first group aligned with the original study while varying `sigma2` and adjusting `lambda0_2` to hold `a0 = 2 * pi * lambda0 * sigma^2` constant (default `a0 = 400` for both groups).

## Structure
- `sim_two_group.R` — helper functions to build scenarios, run simulations, and summarise empirical and approximate RSE for each group.
- Outputs are returned as data frames; redirect them to disk as needed for further plotting (e.g., Figure 2 style).

## Usage
```r
source("EB2Gsim/sim_two_group.R")

plan <- build_two_group_plan(
  sigma1 = 20,
  sigma2 = c(2, 5, 10, 20, 40, 60, 100),
  grids = list(c(6, 6), c(8, 8), c(10, 10)),  # grid sizes; spacing fixed at 2 * sigma1
  a0 = 400,
  D1 = 10,
  D2 = 10,
  nocc = 5,
  detector = "proximity"
)

results <- run_two_group(plan, nrepl = 200, seed = 123, ncores = 4, multisession = TRUE)
head(results$empirical_rse)   # empirical RSE(D-hat) by group and scenario
head(results$approx_counts)   # n, r, and approximation 1/min(sqrt(n), sqrt(r)) by group

# Plot (Figure 2 style)
plot_two_group_rse(results, outfile = "EB2Gsim/fig2_two_group.png")

# Per-replicate outputs with true/estimated parameters and RSE
head(results$per_repl)
```

## Notes
- The trap spacing is always `2 * sigma1`; only grid size changes. This keeps group 1 aligned with the original single-group recommendation.
- When `sigma2` changes, `lambda0_2` is recomputed from `a0` so detectability stays realistic.
- Models are fitted with `D ~ g`, `lambda0 ~ g`, `sigma ~ g` using `groups = "group"` (pattern from `DesignSims1.R`).
- Approximate RSE uses `1 / min(sqrt(n), sqrt(r))`, where `n` is the number of unique animals captured and `r` is the number of recaptures (spatial and non-spatial) per group.
