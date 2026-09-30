# CopMI

CopMI performs Gaussian-copula EM multiple imputation for multivariate
left-censored continuous data. Version: **0.1.0**.

## Installation

R 4.1.0 or later is required. Install the release from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "xhan020322-a11y/CopMI",
  ref = "v0.1.0",
  dependencies = NA,
  upgrade = "never",
  repos = "https://cloud.r-project.org"
)
```

If `remotes` is already installed, omit the first command. Required package
dependencies are installed automatically when missing.

Alternatively, download `CopMI_0.1.0.tar.gz` from the
[v0.1.0 release](https://github.com/xhan020322-a11y/CopMI/releases/tag/v0.1.0)
and install it with dependency resolution:

```r
remotes::install_local(
  "CopMI_0.1.0.tar.gz",
  dependencies = NA,
  upgrade = "never",
  force = TRUE,
  build = FALSE,
  repos = "https://cloud.r-project.org"
)
```

Use the file's full path if it is outside the working directory.
The release archive includes the HTML vignette. A default GitHub source
installation does not rebuild vignettes.

```r
library(CopMI)
```

## Input

Rows are samples and columns are variables. CopMI requires:

- `X_cens`: numeric matrix or data frame containing observed values and either
  the censoring cutoff or `NA` at censored positions;
- `ind`: matching indicator matrix, with `0 = censored` and `1 = observed`;
- `cutoffs`: one censoring cutoff per variable.

Use `input_scale = "log"` when values and cutoffs are already natural-logged.
Use `input_scale = "raw"` for positive original-scale values; completed data
are returned on the same scale supplied by the user.

## Example

The included `nhanes_pah` object contains 1,330 participants and six urinary
PAH variables, with four variables artificially left-censored for illustration.

```r
library(CopMI)
data("nhanes_pah", package = "CopMI")

# Parameter choices:
# input_scale = "log" for already natural-logged data; use "raw" for
# positive original-scale data.
# margin_mode = "normal" when normality after logging is assumed; use
# "select" otherwise to select a marginal family for each variable by BIC.
# m is the number of completed data sets; seed makes the draws reproducible.
# max_iter is the maximum number of EM iterations; tol is the convergence
# tolerance for changes in the copula correlation matrix over three
# consecutive iterations.
# gibbs_burn discards initial sweeps; gibbs_thin separates retained draws.
fit <- copula_em_impute(
  nhanes_pah,
  input_scale = "log",
  margin_mode = "normal",
  m = 5,
  seed = 2026,
  max_iter = 100,
  tol = 1e-4,
  gibbs_burn = 200,
  gibbs_thin = 50
)

# View a concise summary of the fitted imputation model.
summary(fit)

# Extract one completed matrix, all completed matrices, or stacked long data.
completed_1 <- CopMI::complete(fit, action = 1)
completed_all <- CopMI::complete(fit, action = "all")
completed_long <- CopMI::complete(fit, action = "long")
dim(completed_1)
length(completed_all)
head(completed_long)

# Inspect convergence and the fitted copula correlation matrix.
fit_diagnostics <- CopMI::diagnostics(fit)
fit_diagnostics$converged
fit_diagnostics$n_iter
fit_diagnostics$Sigma_hat
fit_diagnostics$sampling$kept_sweeps
```

A single systematic-scan Gibbs chain supplies all completed data sets. Each
sweep updates every censored latent coordinate once. With the settings above,
the five matrices are retained at sweeps 250, 300, 350, 400, and 450.
`diagnostics(fit)$sampling` contains the sweep schedule and two traces: the mean
censored latent value and the mean absolute coordinate update per sweep.
EM convergence and completion of the Gibbs sweeps are separate diagnostics;
neither the sweep count nor these summaries establish adequate mixing.

## Numerical policy

- Marginal candidates use the ordered optimizers on the same censored
  likelihood. Failed candidates are recorded and excluded from BIC selection;
  failure of every candidate stops the fit.
- Positive-definite projection is permitted for the initial pairwise
  correlation estimate only. EM and Gibbs use Cholesky-based conditional
  calculations without ridge, pseudoinverses, or automatic restarts.
- E-step moments come from `tmvtnorm::mtmvnorm`. Invalid moments stop the fit;
  neither untruncated moments nor a Monte Carlo moment-estimation fallback
  is used. The integration library itself may use randomized probability
  integration. Its covariance is symmetrized and checked for positive
  semidefiniteness; relative asymmetry is recorded, not treated as an
  integration-error bound.
- Normal scores, truncated-normal sampling, and inverse transforms use log
  tail probabilities. A failed positive-support quantile is solved from the
  same log-CDF equation on the representable positive range. An unresolved
  quantile, nonfinite result, or material cutoff violation stops imputation.
- Only cutoff equality or a floating-point-sized overshoot may be adjusted
  just below the cutoff. No fixed positive floor is imposed. Adjustments and
  quantile root solves are counted separately in `sampling` and
  `inverse_diagnostics`. Observed entries are never imputed.

`converged` requires three consecutive correlation changes below `tol`; it
does not prove a global likelihood optimum. Gibbs draws condition on fitted
parameters and do not include parameter-estimation uncertainty.

```r
diagnostics(fit)$moment_asymmetry_history
diagnostics(fit)$inverse_diagnostics
diagnostics(fit)$sampling$boundary_adjustments
```

If log-transformed normality is not assumed, use marginal selection:

```r
fit_selected <- copula_em_impute(
  nhanes_pah,
  input_scale = "log",
  margin_mode = "select",
  m = 3,
  seed = 2026,
  max_iter = 100,
  tol = 1e-4
)

# View the selected marginal family for each variable and their BIC values.
selected_diagnostics <- CopMI::diagnostics(fit_selected)
selected_diagnostics$family_selected
selected_diagnostics$bic_table

completed_selected <- CopMI::complete(fit_selected, action = 1)
```

## Support of selected margins

In `margin_mode = "sd_shift"` (or `"select"`), every candidate uses the same
shifted observations and cutoffs. BIC selects the family separately for each
variable. Normal and logistic candidates retain real support: their shifted
imputations may be negative or zero. The other nine candidates require
strictly positive shifted imputations. Both groups require finite censored
values strictly below the corresponding shifted cutoff.

These conditions are checked before undoing the shift. No negative draw from
a normal or logistic margin is replaced merely because of its sign. An invalid
positive-support result is not replaced by a constant. Returned log-analysis
values can be negative for either group; `input_scale = "raw"` returns positive
values after exponentiation. The shift enables positive-support candidates;
it does not assert that every candidate assigns zero probability below the anchor.

The Gaussian latent variables need not be positive. All selected margins
share one copula EM and Gibbs engine, with the same burn-in and retention controls.

## Help

After installation, R help pages describe function arguments and return values;
the vignette is a detailed tutorial, and the package also provides four complete
runnable scripts:

```r
?copula_em_impute                              # Main imputation function
?complete                                      # Extract completed data sets
?diagnostics                                   # Inspect model diagnostics
?nhanes_pah                                    # Example-data structure
vignette("copmi-workflow", package = "CopMI") # Detailed tutorial
system.file("examples", package = "CopMI")    # Location of runnable scripts
```

The package includes separate functions for data validation, marginal fitting,
latent transformation, copula estimation, imputation, completion, and
diagnostics. See the installed help pages for their arguments and return values.

The sampling engine separates truncated-normal draws, preparation of groups
with shared censoring patterns, and advancement of the Gibbs chain. Conditional
covariances must be positive definite. Sampling failures raise errors rather
than modifying covariances or substituting fixed values.

The example data are derived from CDC/NCHS
[NHANES 2015–2016 PAH_I](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PAH_I.htm).
