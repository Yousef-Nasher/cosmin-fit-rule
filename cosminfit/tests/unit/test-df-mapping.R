# cfa_df -------------------------------------------------------------------
check(cfa_df(4, 1) == 2, "cfa_df p=4 m=1")
check(cfa_df(10, 2) == 34, "cfa_df p=10 m=2")
check(cfa_df(30, 1) == 405, "cfa_df p=30 m=1")
check(identical(cfa_df(c(4, 10, 30), c(1, 2, 1)), c(2, 34, 405)), "cfa_df vectorised")
check_warning(cfa_df(3), "df = 0 warns")
check(suppressWarnings(cfa_df(3)) == 0, "df = 0 allowed")
check_error(cfa_df(2), "negative df")
check_error(cfa_df(4.5), "non-integer p")
check_error(cfa_df(0), "p = 0")

# n_moments ----------------------------------------------------------------
check(n_moments(10) == 55, "n_moments lavaan")
check(n_moments(10, "offdiag") == 45, "n_moments offdiag")
check_error(n_moments(10, "other"), "bad convention")

# rmsea_equivalent_table -----------------------------------------------------
tab <- rmsea_equivalent_table()
check(identical(names(tab), c("p", "m", "df", "D", "rmsea_equivalent")), "table columns")
check(nrow(tab) == 8, "table rows")
r4 <- tab[tab$p == 4, ]
check(r4$df == 2 && r4$D == 10, "p=4 df and D")
check(abs(r4$rmsea_equivalent - 0.1789) < 1e-3, "p=4 RMSEA equivalent")
r30 <- tab[tab$p == 30, ]
check(r30$df == 405 && r30$D == 465, "p=30 df and D")
check(abs(r30$rmsea_equivalent - 0.0857) < 1e-3, "p=30 RMSEA equivalent")
check(all(diff(tab$rmsea_equivalent) < 0), "equivalent decreases with p")
tab_off <- rmsea_equivalent_table(p = 4, convention = "offdiag")
check(abs(tab_off$rmsea_equivalent - 0.08 * sqrt(6 / 2)) < 1e-12, "offdiag convention")

# srmr_to_rmsea / rmsea_to_srmr --------------------------------------------
check(abs(srmr_to_rmsea(0.08, 4, 2) - 0.08 * sqrt(5)) < 1e-12, "srmr_to_rmsea value")
s <- seq(0.01, 0.15, by = 0.01)
for (conv in c("lavaan", "offdiag")) {
  back <- rmsea_to_srmr(srmr_to_rmsea(s, 12, cfa_df(12), conv), 12, cfa_df(12), conv)
  check(max(abs(back - s)) < 1e-12, paste("round trip", conv))
  back2 <- srmr_to_rmsea(rmsea_to_srmr(s, 7, 14, conv), 7, 14, conv)
  check(max(abs(back2 - s)) < 1e-12, paste("reverse round trip", conv))
}
check_warning(srmr_to_rmsea(0.05, 3, 0), "df = 0 in mapping warns")
check(is.na(suppressWarnings(srmr_to_rmsea(0.05, 3, 0))), "df = 0 gives NA")
