# CopMI 0.1.0

- Renumber the existing development package (previous internal label 0.4.0)
  as 0.1.0. Algorithms, public function signatures and example data are unchanged.

- Replace the undocumented toy data with `nhanes_pah`: six real NHANES
  2015-2016 exposure variables from the manuscript's 1,330-person sample,
  with a fixed illustrative left-censoring pattern and explicit indicators.
- Ship only censored values, indicators and thresholds, including CSV files.
  Do not ship complete exposure truth or data/censoring generators.
- Remove the evaluation-only `X_true` argument and field. The fourth argument
  of `make_lod_data()` is now `scale`; use named arguments when migrating.
- Update all help, examples, tests and the vignette to the NHANES data.
- Organize the selected teaching workflow into four standalone scripts for
  data preparation, one-call imputation, staged selection, and raw data frames.
- Retain the normal-assumption/selection modes and L-BFGS-B to Nelder-Mead retry.

# Development history (internal label 0.3.0)

* Removed replacement of failed truncated moments with untruncated moments,
  failed samples with starting values, and invalid inverse values with LOD-based
  constants. Such failures now raise contextual errors.
* Positive-support families reject nonpositive observed data; observations are
  never clipped. Shifted normal/logistic families retain their real support.
* Unified all 11 censored marginal likelihoods and optimizers; retained multiple
  starts and L-BFGS-B followed by Nelder-Mead. Positive parameters use log
  coordinates; empirical parameter bounds and probability floors were removed.
* SD-shift uses each column's observed SD, without IQR or pooled-column substitutes.
* Pairwise initialization also retries the same likelihood with two optimizers,
  and errors if both fail. Defaults are min_obs = 2, maxit = 2000, rho_max = 1.
* Only finite initial correlation matrices may be projected to positive
  definiteness; adjustment magnitude is recorded. Removed EM projection,
  eigenvalue clipping, correlation clipping, and conditional covariance ridges.
* Added allow_unconverged = FALSE to imputation APIs. Explicit exploratory
  draws retain their unconverged flag. Summary event counts now sum events.
* Removed the fitdistrplus runtime dependency and duplicated family-fit wrappers.
* Updated help, examples, tests, and the Chinese documentation for these behavior
  changes. Prior results are not expected to remain numerically identical.

# Development history (internal label 0.2.0)

* Introduced staged public APIs, explicit raw/log input scales, the normal-margin
  assumption, heterogeneous marginal selection, and optimizer-attempt diagnostics.
