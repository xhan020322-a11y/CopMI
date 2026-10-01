# CopMI

**Gaussian copula-based multiple imputation for left-censored exposure data.**

CopMI is an R package for imputing multivariate continuous exposures measured
below a limit of detection (LOD). It combines censored marginal distributions
with a Gaussian copula to model dependence between exposures, then generates
multiple completed data sets for subsequent analysis.

- Fit normal margins on the log scale or select among 11 marginal families
  using the Bayesian information criterion (BIC).
- Supply positive original-scale or natural-log-transformed data, with
  completed data returned on the supplied scale and observed values preserved.
- Run the full procedure in one call, or inspect each stage through separate
  fitting functions and diagnostic outputs.

## Method workflow

```text
exposure values + censoring indicators + detection limits
  -> fit censored marginal distributions
  -> transform to latent Gaussian scores
  -> estimate the copula correlation matrix by EM
  -> draw censored values by conditional Gibbs sampling
  -> return multiple completed data sets on the input scale
```

## Installation

Requires **R 4.1.0 or later**. Install the current version from GitHub:

```r
install.packages("remotes")  # Run once if needed
remotes::install_github("xhan020322-a11y/CopMI", ref = "main",
                        upgrade = "never")
```

Required dependencies (`Matrix`, `mvtnorm`, `statmod`, and `tmvtnorm`) are
installed automatically if missing.

## Usage

The included `nhanes_pah` data contain six urinary polycyclic aromatic
hydrocarbon (PAH) metabolites from 1,330 NHANES 2015–2016 participants. Values
are already natural-log-transformed; four variables have artificial left
censoring for illustration.

### Fit the model and extract completed data

```r
library(CopMI)
data("nhanes_pah", package = "CopMI")

fit <- copula_em_impute(
  nhanes_pah,
  input_scale = "log",
  margin_mode = "normal",
  m = 5,
  seed = 2026
)

summary(fit)

# Extract one completed matrix, all five matrices, or stacked long data.
completed_1 <- CopMI::complete(fit, action = 1)
completed_all <- CopMI::complete(fit, action = "all")
completed_long <- CopMI::complete(fit, action = "long")

dim(completed_1)
# [1] 1330    6
length(completed_all)
# [1] 5
```

Each completed matrix retains the original rows, columns, and observed values.
The long-format output adds `.imp` and `.id` to identify the imputation and
original row. Here, `margin_mode = "normal"` assumes normal marginal
distributions after logging.

### Select marginal distributions

Use `margin_mode = "select"` to select a marginal family separately for each
exposure by BIC. Selection is performed on shifted log-scale data; both modes
use the same Gaussian copula model and imputation procedure.

```r
fit_selected <- copula_em_impute(
  nhanes_pah,
  input_scale = "log",
  margin_mode = "select",
  m = 5,
  seed = 2026
)

selected_diagnostics <- CopMI::diagnostics(fit_selected)
selected_diagnostics$family_selected
selected_diagnostics$bic_table
```

### Check model diagnostics

```r
fit_diagnostics <- CopMI::diagnostics(fit)
fit_diagnostics$converged   # EM convergence status
fit_diagnostics$n_iter      # Number of EM iterations
fit_diagnostics$Sigma_hat   # Estimated copula correlation matrix
fit_diagnostics$sampling    # Gibbs sampling schedule and trace summaries
```

The main controls are:

| Argument | Purpose |
| --- | --- |
| `input_scale` | `"log"` for already natural-logged values and cutoffs; `"raw"` for positive original-scale input. |
| `margin_mode` | `"normal"` for normal log-scale margins; `"select"` for BIC-based marginal selection. |
| `m` | Number of completed data sets; default `5`. |
| `seed` | Random seed for reproducible fitting and sampling. |
| `max_iter`, `tol` | EM iteration limit and correlation-change tolerance; defaults `100` and `1e-4`. |
| `gibbs_burn`, `gibbs_thin` | Initial Gibbs sweeps to discard and sweeps between retained draws; defaults `200` and `50`. |

EM convergence and Gibbs mixing should be assessed separately. Imputations
condition on the fitted model parameters; parameter-estimation uncertainty is
not sampled. CopMI returns completed data sets; downstream model fitting and
pooling are performed separately.

## Use your own data

Prepare three objects, with samples in rows and exposure variables in columns:

| Input | Required format |
| --- | --- |
| `X_cens` | Numeric matrix or data frame. Censored cells may contain `NA` or their detection limit. |
| `ind` | Matching indicator matrix: **`0 = censored`, `1 = observed`**. |
| `cutoffs` | Numeric vector with one detection limit per column, in the same order; `NA` is allowed for fully observed columns. |

For positive original-scale concentrations and detection limits:

```r
dat <- make_lod_data(X_cens, ind, cutoffs)

fit_custom <- copula_em_impute(
  dat,
  input_scale = "raw",
  margin_mode = "normal",
  m = 5,
  seed = 2026
)

completed_custom <- CopMI::complete(fit_custom, action = "all")
```

Use `input_scale = "log"` instead when **both values and cutoffs are already
natural-logged**. The censoring indicator determines which cells are imputed.
This input format represents left censoring with one cutoff per variable;
general missingness, right censoring, and observation-specific cutoffs are
not supported.

## Documentation

- [Workflow tutorial](vignettes/copmi-workflow.Rmd): data preparation, model
  fitting, and diagnostics.
- [Runnable examples](inst/examples): data inspection, quick start, staged
  fitting, and raw-data input.
- [Change log](NEWS.md): changes to the package.

After installation, open the function reference in R:

```r
help("copula_em_impute", package = "CopMI")
help("make_lod_data", package = "CopMI")
help("complete", package = "CopMI")
help("diagnostics", package = "CopMI")
```

## Example data

The example data are derived from CDC/NCHS
[NHANES 2015–2016 PAH_I](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PAH_I.htm).
The illustrative censoring thresholds are artificial and do not represent
laboratory detection limits. See `help("nhanes_pah", package = "CopMI")` for
the data structure and source details.

## License

CopMI is released under the [MIT License](LICENSE.md).
