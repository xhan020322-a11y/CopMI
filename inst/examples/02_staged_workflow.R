library(CopMI)
data("nhanes_pah", package = "CopMI")

dat <- make_lod_data(
  X_cens = nhanes_pah$X_cens,
  ind = nhanes_pah$ind,
  cutoffs = nhanes_pah$cutoffs,
  scale = "log"
)

# Step 1: fit and select each margin by BIC after an internal SD shift.
# The default 11 candidates are norm, logis, lnorm, gamma, weibull, exp,
# invgauss, gengamma, llogis, lomax, and burr.
# Output: selected families, parameter estimates, and candidate fit tables.
margins <- copmi_fit_margins(
  dat,
  input_scale = "log",
  margin_mode = "select",
  shift_k = 3,
  optim_methods = c("L-BFGS-B", "Nelder-Mead")
)
print(margins)
print(margins$family_selected)
print(margins$candidate_table)
print(margins$fits[["URXP01"]]$parameters)

# Step 2: transform to latent Gaussian scores and their censoring bounds.
# Censored entries of Z still represent bounds; they are not imputations.
latent <- copmi_transform(margins)
print(latent)
print(head(latent$Z, 3))
print(latent$z_lod)

# Step 3: estimate the copula correlation matrix by EM.
# Output: Sigma_hat, convergence status, iteration count, and matrix changes.
model <- copmi_fit_copula(
  latent, max_iter = 100, tol = 1e-4, seed = 2026
)
print(model)
print(model$Sigma_hat)
print(model$converged)
print(model$n_iter)
print(model$em_change_history)

# Step 4: draw three completed data sets on the original log input scale.
# With seed omitted, sampling continues the random state saved in model.
fit_selected <- copmi_impute(model, m = 3)
print(summary(fit_selected))
completed_selected <- CopMI::complete(fit_selected, action = 1)
print(head(completed_selected, 3))
selected_diagnostics <- CopMI::diagnostics(fit_selected)
print(selected_diagnostics$fallback_records)
