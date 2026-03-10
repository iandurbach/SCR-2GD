# optimalSpacing2G

Two-group extension of the `secrdesign` optimal spacing idea, searching over **absolute trap spacing (m)** to maximise a user-chosen criterion based on expected captures/recaptures for two groups with their own `D`, `lambda0`, `sigma`, and detection function.

## What it does
- Scales a baseline trap layout (`traps0`) to candidate spacings (metres), rebuilds a mask, and for each group computes expected captures (`En`) and recaptures (`Er`) with `Enrm`.
- Objective options:
  - `criterion = "sum_min"`: maximise `min(En1 + En2, Er1 + Er2)`.
  - `criterion = "all_min"`: maximise `min(En1, En2, Er1, Er2)`.
- Poisson only; no simulation block. Uses `secrdesign:::dfcast`/`replacedefaults` internally to harmonise detection functions/pars.
- Returns a data frame of per-spacing values and the spacing that maximises the chosen criterion.

## Usage
```r
library(secrdesign)
source("optimalSpacing2G/optimalSpacing2G.R")

grid <- make.grid(7, 7, spacing = 50)  # baseline traps (multi-catch)

# One-group reference (secrdesign native)
os1 <- optimalSpacing(
  D = 5,
  traps = grid,
  detectpar = list(lambda0 = 0.2, sigma = 40),
  noccasions = 5,
  plt = FALSE
)
os1$rotRSE$optimum.spacing  # e.g., ~57 m

# Two groups; identical params should recover the same spacing
os2 <- optimalSpacing2G(
  D1 = 5, D2 = 5,
  traps0 = grid,
  detectpar1 = list(lambda0 = 0.2, sigma = 40),
  detectpar2 = list(lambda0 = 0.2, sigma = 40),
  noccasions = 5,
  spacing_m = seq(40, 120, 1),
  criterion = "sum_min"
)
os2$optimum.spacing         # ~57 m, matches one-group result
head(os2$values)            # columns: spacing, En1, En2, Er1, Er2, crit
```

## Notes
- `spacing_m` sets the candidate spacings (metres); choose the range/resolution you need. The optimiser interpolates; if you want the best evaluated spacing, replace `interpCritMax` with a `which.max` on `values$crit`.
- `mask` buffer defaults to `xsigma * max(sigma1, sigma2) + spacing`; adjust inside `getCrit2G_abs` if you prefer a different rule.
- Trap scaling is done by multiplying `x`/`y` and spacing attributes by the requested factor; detector type is preserved.***
