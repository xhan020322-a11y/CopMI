library(CopMI)

# Read the already-censored concentrations (ng/L), indicators, and thresholds.
# Censored CSV entries are NA. No complete exposure truth is included.
pah_table <- read.csv(
  system.file("extdata", "nhanes_pah_censored.csv", package = "CopMI"),
  row.names = "SEQN"
)
metadata <- read.csv(
  system.file("extdata", "nhanes_pah_metadata.csv", package = "CopMI")
)

# Both data frames must have matching row/column names and variable order.
# ind_raw: 0 = left-censored, 1 = observed; lod_raw: one raw cutoff per column.
X_raw <- pah_table[metadata$variable]
ind_raw <- pah_table[paste0(metadata$variable, "_ind")]
names(ind_raw) <- metadata$variable
lod_raw <- metadata$cutoff_raw
print(head(X_raw, 3))
print(head(ind_raw, 3))
print(lod_raw)

# Pass a numeric data.frame directly, together with indicators and cutoffs.
# "raw" logs values internally; "normal" assumes normality after logging.
fit_raw <- copula_em_impute(
  x = X_raw,
  ind = ind_raw,
  cutoffs = lod_raw,
  input_scale = "raw",
  margin_mode = "normal",
  optim_methods = c("L-BFGS-B", "Nelder-Mead"),
  m = 3,
  seed = 2026
)

# Output is a completed matrix on the supplied raw scale (ng/L).
completed_raw <- CopMI::complete(fit_raw, action = 1)
print(summary(fit_raw))
print(dim(completed_raw))
print(head(completed_raw, 3))
