# CopMI

CopMI performs Gaussian-copula EM multiple imputation for multivariate
left-censored continuous data.

## Features

- Normal margins or BIC-based selection among 11 marginal families.
- Raw or log-scale input, with completed data returned on the supplied scale.
- A complete imputation workflow, separate fitting functions, and model diagnostics.

## Installation

Requires R 4.1.0 or later. Missing dependencies are installed automatically.

```r
install.packages("remotes")  # Run once if needed
remotes::install_github("xhan020322-a11y/CopMI", ref = "v0.1.0",
                        upgrade = "never")
```

A source installation archive is also available on the
[release page](https://github.com/xhan020322-a11y/CopMI/releases/tag/v0.1.0).

## Quick example

The included `nhanes_pah` data provide an example with illustrative left censoring.

```r
library(CopMI)
data("nhanes_pah", package = "CopMI")

fit <- copula_em_impute(nhanes_pah, margin_mode = "normal", m = 5, seed = 2026)
completed <- CopMI::complete(fit, "all")
summary(fit)
```

This example assumes normal margins on the log scale. Use `margin_mode = "select"`
for BIC-based marginal selection.

For your own data, use `make_lod_data()` to supply the values, censoring indicators
(`0 = censored`, `1 = observed`), and one cutoff per variable.

## Documentation

- [Workflow tutorial](vignettes/copmi-workflow.Rmd)
- [Runnable examples](inst/examples)
- R help: `?copula_em_impute`, `?make_lod_data`, `?complete`, `?diagnostics`

Example data are derived from CDC/NCHS
[NHANES 2015–2016 PAH_I](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PAH_I.htm).
