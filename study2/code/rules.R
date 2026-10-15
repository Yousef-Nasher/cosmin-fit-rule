# rules.R -- shared by the Study 2 analysis scripts: which indices each analysis arm
# gives a reader, and the rating rules applied to them.

# Clean fits: converged, with admissible estimates.
clean_fits <- function(D) D[D$status == "ok" & D$admiss %in% TRUE, ]

# Indices by arm. For WLSMV the scaled indices are primary and the robust indices are the
# sensitivity analysis; the ML arms have one set of indices.
idx <- function(E, version) {
  w <- E$arm == "wlsmv_ord"
  pick <- function(ml, sc, rob) ifelse(w, if (version == "robust") rob else sc, ml)
  data.frame(cfi = pick(E$cfi, E$cfi.sc, E$cfi.rob), tli = pick(E$tli, E$tli.sc, E$tli.rob),
             rmsea = pick(E$rmsea, E$rmsea.sc, E$rmsea.rob), srmr = E$srmr)
}

# A missing index counts as not passing.
f <- function(v) !is.na(v) & v

rules <- function(X, df, Dd) {
  c_ <- f(X$cfi > .95); t_ <- f(X$tli > .95); r_ <- f(X$rmsea < .06); s_ <- f(X$srmr < .08)
  s_sc <- f(X$srmr < .06 * sqrt(df / Dd))            # SRMR cut-off scaled by sqrt(df / D)
  s09 <- f(X$srmr <= .09); r06 <- f(X$rmsea <= .06)
  data.frame(
    COSMIN          = c_ | t_ | r_ | s_,
    COSMIN_noRMSEA  = c_ | t_ | s_,
    SRMR_df_scaled  = c_ | t_ | r_ | s_sc,
    two_of_three    = (c_ + r_ + s_) >= 2,
    all_three       = c_ & r_ & s_,
    HuBentler       = (f(X$cfi >= .95) & s09) | (r06 & s09),   # two of the Hu-Bentler combination rules
    HuBentler_CFI96 = (f(X$cfi >= .96) & s09) | (r06 & s09),
    HuBentler_TLI96 = (f(X$tli >= .96) & s09) | (r06 & s09),
    CFI_only = c_, RMSEA_only = r_, SRMR_only = s_,
    rmsea_rescue = r_ & !c_ & !t_ & !s_,               # RMSEA is the only clause that passes
    srmr_sole    = s_ & !c_ & !t_ & !r_,
    cfitli_sole  = (c_ | t_) & !r_ & !s_,
    srmr_gt_rmsea = f(X$srmr > X$rmsea))
}

mcse <- function(x) sqrt(mean(x) * (1 - mean(x)) / length(x))   # Monte Carlo standard error

# Read the unzipped simulation output and attach the population severity flags.
read_study2 <- function() {
  fs <- list.files("out/final", pattern = "\\.csv$", full.names = TRUE)
  D  <- do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
  P  <- read.csv("out/population.csv", stringsAsFactors = FALSE)
  PU <- P[P$est == "ULS", ]; PM <- P[P$est == "ML", ]
  D$cond <- sub("_n[0-9]{4}$", "", D$cell)
  D$conseq     <- PU$consequential[match(D$cond, PU$cond)]   # primary severity flag (ULS)
  D$conseq_ml  <- PM$consequential[match(D$cond, PM$cond)]   # sensitivity
  D$omega_over <- PU$omega_over[match(D$cond, PU$cond)]
  D$Dd <- D$p * (D$p + 1) / 2
  attr(D, "n_files") <- length(fs)
  D
}
