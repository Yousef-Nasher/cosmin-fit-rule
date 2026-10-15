# table_rules.R -- Study 2: false-fail and false-pass rates of every rule, by analysis arm
# and N, with Monte Carlo standard errors (source of the rules table and of Figure 1).
# Run from study2/ after analysis.R:  Rscript code/table_rules.R
source("code/rules.R")
X <- readRDS("out/study2_analysis.rds"); E <- X$E
out <- list()
for (ver in c("scaled", "robust")) {
  R <- rules(idx(E, ver), E$df, E$Dd)
  for (arm in unique(E$arm)) for (n in sort(unique(E$n))) {
    if (ver == "robust" && arm != "wlsmv_ord") next
    s <- E$arm == arm & E$n == n
    cor <- s & E$kind == "none"; con <- s & E$conseq
    for (r in c("COSMIN", "COSMIN_noRMSEA", "SRMR_df_scaled", "two_of_three", "all_three", "HuBentler"))
      out[[length(out) + 1]] <- data.frame(version = ver, arm = arm, n = n, rule = r,
        false_fail = 100 * mean(!R[[r]][cor]), mcse_ff = 100 * sqrt(mean(!R[[r]][cor]) * mean(R[[r]][cor]) / sum(cor)),
        false_pass = 100 * mean(R[[r]][con]), mcse_fp = 100 * sqrt(mean(R[[r]][con]) * mean(!R[[r]][con]) / sum(con)),
        n_correct = sum(cor), n_conseq = sum(con))
  }
}
T <- do.call(rbind, out)
write.csv(T, "out/table_rules_by_N.csv", row.names = FALSE)
print(T[T$version == "scaled" & T$arm == "wlsmv_ord", ], digits = 3, row.names = FALSE)
