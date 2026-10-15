# rate_cfa must agree exactly with the reference used in the simulation study.
ref_rules <- function(X, df, Dd) {
  c_ <- X$cfi > .95; t_ <- X$tli > .95; r_ <- X$rmsea < .06; s_ <- X$srmr < .08
  s_sc <- X$srmr < .06 * sqrt(df / Dd)
  data.frame(COSMIN = c_ | t_ | r_ | s_, COSMIN_noRMSEA = c_ | t_ | s_,
    SRMR_df_scaled = c_ | t_ | r_ | s_sc, two_of_three = (c_ + r_ + s_) >= 2,
    all_three = c_ & r_ & s_,
    HuBentler = (X$cfi >= .95 & X$srmr <= .09) | (X$rmsea <= .06 & X$srmr <= .09))
}

set.seed(1)
n <- 5000
X <- data.frame(cfi = runif(n, .85, 1), tli = runif(n, .85, 1),
                rmsea = runif(n, 0, .15), srmr = runif(n, 0, .15))
X$p <- sample(4:20, n, replace = TRUE)
X$df <- cfa_df(X$p)
ref <- ref_rules(X, X$df, n_moments(X$p))
map <- c(COSMIN = "cosmin", COSMIN_noRMSEA = "no_rmsea",
         SRMR_df_scaled = "srmr_df", two_of_three = "two_of_three",
         all_three = "all_three", HuBentler = "hu_bentler")

for (nm in names(map)) {
  got <- rate_cfa(X$cfi, X$tli, X$rmsea, X$srmr, rule = map[[nm]],
                  p = X$p, df = X$df)
  want <- ifelse(ref[[nm]], "sufficient", "insufficient")
  check(identical(got, want), paste("matches reference:", nm))
}
check(var(ref$COSMIN) > 0 && var(ref$all_three) > 0, "reference has both outcomes")

# With complete data, "reported" and "unknown" agree for every rule.
all_rules <- cosmin_rules()$rule
a <- rerate(X, rules = all_rules, missing = "reported")
b <- rerate(X, rules = all_rules, missing = "unknown")
for (r in all_rules) {
  col <- paste0("rating_", r)
  check(identical(a[[col]], b[[col]]), paste("reported == unknown:", r))
  check(!any(a[[col]] == "indeterminate"), paste("no indeterminate:", r))
}
