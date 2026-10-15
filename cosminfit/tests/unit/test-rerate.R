d <- data.frame(cfi = c(0.97, 0.91, NA, NA, 0.90),
                tli = c(NA, 0.89, NA, NA, NA),
                rmsea = c(0.05, 0.05, NA, 0.08, 0.07),
                srmr = c(0.04, 0.09, 0.07, NA, 0.05),
                p = c(10, 6, 8, 5, 6), df = c(35, 9, 20, 5, 9))
out <- rerate(d)
check(all(names(d) %in% names(out)), "original columns kept")
check(all(paste0("rating_", c("cosmin", "no_rmsea", "srmr_df", "two_of_three",
                              "all_three")) %in% names(out)), "rating columns")
check(identical(out$n_indices, c(3, 4, 1, 1, 3)), "n_indices")
check(identical(out$rmsea_rescue, c(FALSE, TRUE, NA, FALSE, FALSE)), "rmsea_rescue")
check(identical(out$rating_cosmin,
                c("sufficient", "sufficient", "sufficient", "insufficient", "sufficient")),
      "cosmin column")
check(identical(out$rating_srmr_df[3], "insufficient"), "srmr_df row 3 (.07 > .0447)")
check(identical(out$rating_all_three,
                c("sufficient", "insufficient", "indeterminate", "indeterminate", "insufficient")),
      "all_three column")
for (r in c("cosmin", "all_three")) {
  check(identical(out[[paste0("rating_", r)]],
                  rate_cfa(d$cfi, d$tli, d$rmsea, d$srmr, rule = r, p = d$p, df = d$df)),
        paste("rerate equals rate_cfa:", r))
}

# Custom and partial column names; absent columns count as not reported.
d2 <- data.frame(CFI = c(0.97, 0.90), SRMR = c(0.10, 0.10))
out2 <- rerate(d2, rules = "cosmin", cols = list(cfi = "CFI", srmr = "SRMR"))
check(identical(out2$rating_cosmin, c("sufficient", "insufficient")), "custom cols")
check(identical(out2$n_indices, c(2, 2)), "custom cols n_indices")
check(all(is.na(out2$rmsea_rescue)), "no RMSEA column -> rescue NA")
out3 <- rerate(d2, rules = "cosmin", missing = "unknown",
               cols = list(cfi = "CFI", srmr = "SRMR"))
check(identical(out3$rating_cosmin, c("sufficient", "indeterminate")), "unknown mode")
check_error(rerate(d, rules = "nope"), "bad rule")
check_error(rerate(as.list(d)), "data must be data.frame")
