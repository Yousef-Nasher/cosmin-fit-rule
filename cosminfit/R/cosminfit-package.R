#' cosminfit: rate CFA fit against the COSMIN criterion
#'
#' Rates confirmatory factor analysis (CFA) results against the COSMIN
#' structural validity criterion "CFI or TLI or comparable measure > 0.95
#' OR RMSEA < 0.06 OR SRMR < 0.08" and against alternative rules, and gives
#' the population SRMR-RMSEA mapping.
#'
#' Main functions: [rate_cfa()], [rerate()], [cosmin_rules()],
#' [srmr_to_rmsea()], [rmsea_equivalent_table()], [cfa_df()],
#' [wilson_ci()], [cluster_boot_prop()].
#'
#' @references
#' Prinsen CAC, Mokkink LB, Bouter LM, et al. (2018). COSMIN guideline for
#' systematic reviews of patient-reported outcome measures. Quality of Life
#' Research, 27, 1147-1157.
#'
#' Prinsen CAC, Vohra S, Rose MR, et al. (2016). Guideline for selecting
#' outcome measurement instruments for outcomes included in a Core Outcome
#' Set. Trials, 17, 449.
#'
#' Mokkink LB, Elsman EBM, Terwee CB (2024). COSMIN guideline for systematic
#' reviews of patient-reported outcome measures version 2.0. Quality of Life
#' Research, 33, 2929-2939.
#'
#' Shi D, Maydeu-Olivares A, Rosseel Y (2020). Assessing fit in ordinal
#' factor analysis models: SRMR vs. RMSEA. Structural Equation Modeling, 27,
#' 1-15.
#' @keywords internal
"_PACKAGE"
