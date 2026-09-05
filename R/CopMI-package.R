#' CopMI: Gaussian copula imputation for left-censored exposures
#'
#' Use [copula_em_impute()] for the complete workflow or call
#' [make_lod_data()], [copmi_fit_margins()], [copmi_transform()],
#' [copmi_fit_copula()], and [copmi_impute()] separately. Use [complete()] to
#' extract data and [diagnostics()] to inspect convergence and optimizer retries.
#' @section Modelling choices:
#' `input_scale` states whether exposures are raw positive values or already
#' logged. `margin_mode = "normal"` assumes normal margins on the log scale;
#' `"sd_shift"` (alias `"select"`) uses BIC selection after an SD shift.
#' The default optimizer order is L-BFGS-B then Nelder-Mead on log-transformed positive parameters.
#' @section Scope:
#' The implementation conditions imputations on estimated margins and copula
#' correlation. It does not sample model-parameter uncertainty or pool downstream
#' regression analyses. It supports left censoring with one cutoff per column.
#' @seealso [nhanes_pah], `vignette("copmi-workflow", package = "CopMI")`
#' @keywords internal
"_PACKAGE"
