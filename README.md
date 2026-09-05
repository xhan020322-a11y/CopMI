# CopMI

CopMI performs Gaussian-copula EM multiple imputation for multivariate
left-censored continuous data. Version: **0.1.0**.

## Installation

```r
install.packages("remotes")
remotes::install_github("xhan020322-a11y/CopMI")
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
# tolerance for changes in the copula correlation matrix.
fit <- copula_em_impute(
  nhanes_pah,
  input_scale = "log",
  margin_mode = "normal",
  m = 3,
  seed = 2026,
  max_iter = 100,
  tol = 1e-4
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

The example data are derived from CDC/NCHS
[NHANES 2015–2016 PAH_I](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PAH_I.htm).
