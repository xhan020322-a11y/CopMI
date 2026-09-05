library(CopMI)
data("nhanes_pah", package = "CopMI")

dat <- make_lod_data(
  X_cens = nhanes_pah$X_cens,
  ind = nhanes_pah$ind,
  cutoffs = nhanes_pah$cutoffs,
  scale = "log"
)

# Assume normal margins after logging; this is a modelling choice.
# Nelder-Mead retries the same likelihood only if L-BFGS-B fails.
# Return a copmi_mi object with three completed data sets and diagnostics.
fit <- copula_em_impute(
  dat,
  input_scale = "log",
  margin_mode = "normal",
  optim_methods = c("L-BFGS-B", "Nelder-Mead"),
  m = 3,
  seed = 2026,
  max_iter = 100,
  tol = 1e-4
)
print(fit)
print(summary(fit))

# Extract one matrix, all three matrices, their cellwise mean, or stacked rows.
# Outputs retain the log input scale. The mean is descriptive, not MI pooling.
completed_log <- CopMI::complete(fit, action = 1)
completed_all <- CopMI::complete(fit, action = "all")
completed_mean <- CopMI::complete(fit, action = "mean")
completed_long <- CopMI::complete(fit, action = "long")
print(head(completed_log, 3))
print(length(completed_all))
print(dim(completed_mean))
print(head(completed_long, 3))

# Inspect EM convergence, fitted families, optimization attempts, and retries.
fit_diagnostics <- CopMI::diagnostics(fit)
print(names(fit_diagnostics))
print(fit_diagnostics$converged)
print(fit_diagnostics$n_iter)
print(fit_diagnostics$family_selected)
print(head(fit_diagnostics$optimization))
print(fit_diagnostics$fallback_records)
