# hu_bentler_sensitivity.R -- Study 2: sensitivity of the Hu-Bentler combination rule to
# its CFI/TLI cut-off.
# Run from study2/:  Rscript code/hu_bentler_sensitivity.R > results/hu_bentler_sensitivity.txt
# Hu and Bentler (1999) proposed two-index rules in which a model is rejected if either
# index misses its cut-off. The rules below join two of those combinations with "or":
#   HB_CFI95 : (CFI >= .95 & SRMR <= .09) | (RMSEA <= .06 & SRMR <= .09)   (the rule of the main analysis)
#   HB_CFI96 : (CFI >= .96 & SRMR <= .09) | (RMSEA <= .06 & SRMR <= .09)
#   HB_TLI96 : (TLI >= .96 & SRMR <= .09) | (RMSEA <= .06 & SRMR <= .09)
options(width = 150)
source("code/rules.R")
D <- read_study2(); E <- clean_fits(D)
cat(sprintf("job files %d | rows %d | clean fits %d\n", attr(D, "n_files"), nrow(D), nrow(E)))

# false-fail = % of correct models failing the rule
# false-pass = % of consequentially misspecified models passing the rule
table_rates <- function(E, version, arms) {
  R <- rules(idx(E, version), E$df, E$Dd)[, c("COSMIN", "HuBentler", "HuBentler_CFI96", "HuBentler_TLI96")]
  names(R) <- c("COSMIN", "HB_CFI95", "HB_CFI96", "HB_TLI96")
  out <- NULL
  for (a in arms) for (n in sort(unique(E$n))) {
    g0 <- E$arm == a & E$n == n & E$kind == "none"
    g1 <- E$arm == a & E$n == n & E$conseq
    row <- data.frame(arm = a, N = n, n_correct = sum(g0), n_conseq = sum(g1))
    for (r in names(R)) {
      row[[paste0(r, "_ff")]] <- round(100 * mean(!R[[r]][g0]), 1)
      row[[paste0(r, "_fp")]] <- round(100 * mean(R[[r]][g1]), 1)
    }
    out <- rbind(out, row)
  }
  out
}

arms <- c("ml_cont", "ml_ord_as_cont", "wlsmv_ord")
cat("\nff = false-fail % (correct models); fp = false-pass % (consequential misspecification)\n")
cat("HB_CFI95: (CFI>=.95 & SRMR<=.09) | (RMSEA<=.06 & SRMR<=.09)\n")
cat("HB_CFI96: (CFI>=.96 & SRMR<=.09) | (RMSEA<=.06 & SRMR<=.09)\n")
cat("HB_TLI96: (TLI>=.96 & SRMR<=.09) | (RMSEA<=.06 & SRMR<=.09)\n")
cat("The COSMIN columns repeat the main analysis for comparison.\n")

cat("\n=== PRIMARY: WLSMV scaled indices and both ML arms ===\n")
print(table_rates(E, "scaled", arms), row.names = FALSE)

cat("\n=== SENSITIVITY: WLSMV robust indices ===\n")
print(table_rates(E, "robust", "wlsmv_ord"), row.names = FALSE)
