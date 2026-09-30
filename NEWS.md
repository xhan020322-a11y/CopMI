# CopMI 0.1.0

- Gaussian-copula EM multiple imputation for multivariate left-censored data.
- Raw and log input scales, normal margins, and BIC selection among 11
  marginal families after an SD-based shift.
- Persistent systematic-scan Gibbs sampling with explicit burn-in and
  retention controls; defaults are 200 and 50 sweeps, respectively.
- EM convergence requiring three consecutive correlation changes below
  the tolerance, with conditional matrices reused by censoring pattern.
- Log-tail sampling, inverse-CDF validation, and diagnostics for covariance,
  marginal support, and floating-point boundary adjustments.
- Staged fitting interfaces, NHANES example data, runnable examples,
  function documentation, and a workflow vignette.
