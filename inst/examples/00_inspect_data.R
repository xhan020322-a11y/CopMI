library(CopMI)
data("nhanes_pah", package = "CopMI")

# 1,330 participants and six PAH variables; four variables are left-censored.
# Values and cutoffs are already natural-logged. ind: 0 = censored, 1 = observed.
X <- nhanes_pah$X_cens
ind <- nhanes_pah$ind
lod <- nhanes_pah$cutoffs
print(dim(X))
print(head(X, 3))
print(head(ind, 3))

censor_summary <- data.frame(
  variable = colnames(X),
  n_censored = colSums(ind == 0L),
  percent_censored = round(100 * colMeans(ind == 0L), 2),
  cutoff_log = lod,
  row.names = NULL
)
print(censor_summary)

# Validate the three inputs and return a copmi_lod_data object.
# Censored entries represent a bound, not a measured concentration.
dat <- make_lod_data(X_cens = X, ind = ind, cutoffs = lod, scale = "log")
print(dat)
