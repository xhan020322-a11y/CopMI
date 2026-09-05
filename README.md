# CopMI

Gaussian copula EM imputation for multivariate left-censored exposure data.
Version 0.1.0 separates validation, marginal fitting, latent transformation,
copula estimation, conditional draws, and diagnostic extraction.

## Install

Install the runtime dependencies once, then install the supplied source archive:

```r
install.packages(c("Matrix", "mvtnorm", "statmod", "tmvtnorm"))
install.packages("CopMI_0.1.0.tar.gz", repos = NULL, type = "source")
library(CopMI)
```

This is a local source package; it has not been submitted to CRAN.
From an unpacked source directory, use `R CMD INSTALL CopMI`.

## Choose the input scale and modelling assumption

| Input data | Assumption | Parameters |
|---|---|---|
| Already natural-logged values and cutoffs | Normal log margins | `input_scale = "log", margin_mode = "normal"` |
| Positive original values and cutoffs | Normal after natural log | `input_scale = "raw", margin_mode = "normal"` |
| Already logged values and cutoffs | Select a family per column | `input_scale = "log", margin_mode = "select"` |
| Positive original values and cutoffs | Select on shifted log scale | `input_scale = "raw", margin_mode = "select"` |

`"select"` is an alias for the original default `"sd_shift"`. There is no
automatic normality test. Completed matrices use the scale supplied by the
caller: raw in, raw out; logged in, logged out. The `scale` label in the data
object is descriptive; the explicit `input_scale` argument controls conversion.

Each marginal fit first uses **L-BFGS-B**. If it errors, fails to converge,
or returns invalid results, it retries the **same censored likelihood and
log-parameterization with Nelder-Mead**. The default is
`optim_methods = c("L-BFGS-B", "Nelder-Mead")` in both modes.

## NHANES example data

`nhanes_pah` contains six urinary PAH metabolites from the manuscript's 1,330
NHANES 2015-2016 participants. The concentrations are real and the censoring is
artificial, fixed for this teaching example. URXP01 is censored below 346.9 ng/L
(133 values, 10%); URXP03 below 43 ng/L (396 values, 29.77%); URXP04 below
80 ng/L (265 values, 19.92%); and URXP06 below 87 ng/L (529 values, 39.77%).
URXP02 and URXP25 remain observed. These are illustrative thresholds, not laboratory
LODs or a particular evaluation replicate. This example does not establish
normality after logging.

The packaged R object is already on the natural-log scale. Its `ind` matrix uses
**0 = censored, 1 = observed**; `X_cens` stores log(LOD) placeholders in censored
cells. Raw-scale CSV data with NA in censored cells and matching `_ind` columns
are installed under `system.file("extdata", package = "CopMI")`, together with a
variable/threshold table. See `?nhanes_pah`. Neither complete exposure truth nor
data-generation/censoring-generation code is shipped.

Source: CDC/NCHS, [NHANES 2015-2016 PAH_I](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/PAH_I.htm),
linked to DEMO_I, DPQ_I, BMX_I and ALB_CR_I for cohort selection. Original public
files are freely available from NHANES. This derived example is not endorsed by
CDC, HHS or the US Government. Examples illustrate imputation and do not perform
survey-weighted population inference.

## One-call workflow

```r
library(CopMI)
data("nhanes_pah", package = "CopMI")

dat <- make_lod_data(
  X_cens = nhanes_pah$X_cens,
  ind = nhanes_pah$ind,
  cutoffs = nhanes_pah$cutoffs,
  scale = "log"
)

# The example values are used as already logged values.
fit <- copula_em_impute(
  dat,
  input_scale = "log",
  margin_mode = "normal",
  optim_methods = c("L-BFGS-B", "Nelder-Mead"),
  m = 3,
  max_iter = 100,
  seed = 2026
)

summary(fit)
completed_log <- CopMI::complete(fit, action = 1)
completed_all <- CopMI::complete(fit, action = "all")
completed_mean <- CopMI::complete(fit, action = "mean")
completed_long <- CopMI::complete(fit, action = "long")
fit_diagnostics <- CopMI::diagnostics(fit)
fit_diagnostics$converged
fit_diagnostics$family_selected
head(fit_diagnostics$optimization)
fit_diagnostics$fallback_records
```

For heterogeneous selection, set `margin_mode = "select"`. The default candidate
set is `norm`, `logis`, `lnorm`, `gamma`, `weibull`, `exp`, `invgauss`,
`gengamma`, `llogis`, `lomax`, and `burr`. Use `margin_candidates` to narrow it.

## Separate steps with explicit outputs

```r
# Continue with dat above; select margins instead of assuming normality.
margins <- copmi_fit_margins(
  dat, input_scale = "log", margin_mode = "select", shift_k = 3,
  optim_methods = c("L-BFGS-B", "Nelder-Mead")
)
margins$family_selected
margins$candidate_table
latent <- copmi_transform(margins)
model <- copmi_fit_copula(latent, max_iter = 100, tol = 1e-4, seed = 2026)
model$Sigma_hat
fit_selected <- copmi_impute(model, m = 3)
completed_selected <- CopMI::complete(fit_selected, action = 1)
```

| Function | Responsibility | Main output |
|---|---|---|
| `make_lod_data()` | Validate values, indicators, cutoffs and names | `copmi_lod_data`: `X_cens`, `ind`, `cutoffs` |
| `copmi_shift()` | Compute and apply an SD shift on a supplied analysis scale | `copmi_shift`: `dat`, `anchor`, `sd` |
| `copmi_fit_margin()` | Fit/compare a single column on the supplied scale | `copmi_margin`: `parameters`, `candidate_table`, `optimization` |
| `copmi_fit_margins()` | Prepare scale and fit every column | `copmi_margins`: `fits`, `family_selected`, tables |
| `copmi_transform()` | Convert fitted marginal probabilities to Gaussian scores | `copmi_latent`: `Z`, `ind`, `z_lod` |
| `copmi_fit_copula()` | Initialize and estimate correlation by EM | `copmi_copula`: `Sigma_hat`, convergence and history |
| `copmi_impute()` | Draw and back-transform completed data | `copmi_mi`: `imp_list`, diagnostics |
| `copula_em_impute()` | Coordinate the full workflow | `copmi_mi` |
| `complete()` | Extract a matrix, list, mean, or long data frame | Specified data representation |
| `diagnostics()` | Extract fitting/retry/EM/sampling information | Named list of tables and values |
| `print()`, `summary()` | Display/describe supported S3 objects | Original object invisibly / one-row summary |

Every public function has installed R help with arguments, return fields, and
executable examples: `?copmi_fit_margin`, `?copula_em_impute`, etc.
Read `vignette("copmi-workflow", package = "CopMI")` for a complete tutorial.
Runnable scripts are installed under `system.file("examples", package = "CopMI")`.

## Standalone examples

Each script starts with `library(CopMI)`, loads its own inputs, and leaves the
results in your workspace. All use the full censored NHANES sample; imputation
examples return three completed data sets with `seed = 2026`.

| Script | What it shows |
|---|---|
| [00_inspect_data.R](inst/examples/00_inspect_data.R) | Values, censor indicators, thresholds, a censoring summary, and `make_lod_data()` |
| [01_quick_start.R](inst/examples/01_quick_start.R) | Normal-assumption one-call imputation, four `complete()` formats, and `diagnostics()` |
| [02_staged_workflow.R](inst/examples/02_staged_workflow.R) | Select among 11 marginal families, transform, fit the copula, and impute; inspect every stage |
| [03_raw_data_frame.R](inst/examples/03_raw_data_frame.R) | Read the censored CSV and pass raw-scale data frames, indicators, and cutoffs directly |

For example, run the complete one-call script after installation:

```r
source(system.file("examples", "01_quick_start.R", package = "CopMI"))
```

To read a script before running it, use
`file.show(system.file("examples", "01_quick_start.R", package = "CopMI"))`.

## Data and inference scope

- Rows are samples, columns are exposures. Indicators are exactly `1` observed
  and `0` left-censored. Censored placeholders can be `NA` and are standardized.
- Values and cutoffs must share the supplied scale. Each censored column needs
  a finite cutoff. An observed value cannot fall below its cutoff.
- Marginal fitting requires at least two observed values per column;
  non-exponential families additionally require observed variation. Fully censored columns are not estimable by this interface.
- The SD-shift mode uses `anchor = cutoff - shift_k * observed_sd` (or the
  minimum observed value in place of an absent cutoff). Positive-support families imply values above the anchor. Normal/logistic
  families retain real support and can generate values below the anchor.
- Imputations condition on fitted marginal parameters and a fitted correlation
  matrix. Parameter and family-selection uncertainty are not sampled. The
  Gibbs sampler starts inside the truncation region; assess mixing
  and inferential performance for the intended application.
- `complete(fit, "mean")` is descriptive. It is not a substitute for analyzing
  the completed data sets separately and accounting for imputation uncertainty.
- `converged` describes the EM matrix-change criterion only. Inspect the
  candidate, optimizer, initialization, and fallback tables separately.

## Failure policy in 0.1.0

Positive parameters are optimized on the log scale, without empirical upper
bounds. CDFs and densities use log tails; probabilities and observations are
not clipped. Failed E-steps, draws, or inverse transforms raise contextual
errors instead of substituting untruncated moments, starting values, or LOD/sqrt(2).
Only the initial pairwise correlation matrix may be projected to positive
definiteness, with its adjustment recorded. EM updates do not add ridges.

Imputation requires a converged copula by default. For an explicitly exploratory
run, `allow_unconverged = TRUE` returns draws with `converged = FALSE` preserved.
SD-shift uses only each variable's observed SD. Old numerical results are not
promised to match: the former substitutions and parameter bounds changed them.

## Development

```r
install.packages(c("roxygen2", "testthat", "knitr", "rmarkdown"))
roxygen2::roxygenise("CopMI")
testthat::test_local("CopMI")
```

```sh
R CMD build CopMI
R CMD check --as-cran CopMI_0.1.0.tar.gz
```

Building the HTML vignette needs Pandoc. A PDF reference manual additionally
needs a working LaTeX installation; `--no-manual` skips only that PDF check.
Package layout, namespace registration, Rd help, examples, tests, and vignettes
follow the [Writing R Extensions manual](https://cran.r-project.org/doc/manuals/r-release/R-exts.html).
