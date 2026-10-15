# analysis2b.R -- Study 2b: analysis as fixed in plan.md.
# Run from study2b/:  Rscript code/analysis2b.R > results/study2b_results.txt
# Reads the unzipped simulation output (out/chunks/*.csv) and out/population2b.csv.
options(width = 180)
fs <- list.files("out/chunks", "\\.csv$", full.names = TRUE)
D  <- do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
P  <- read.csv("out/population2b.csv"); PU <- P[P$est == "ULS", ]
D$conseq <- PU$consequential[match(paste(D$p, D$load, D$kind), paste(PU$p, PU$load, PU$kind))]
D$fmt <- factor(D$ct, levels = c("c2_sym","c2_skew","c3_sym","c3_skew","c5_sym","c5_skew"))
reps <- tapply(D$rep, D$cell, function(x) length(unique(x)))
cat(sprintf("jobs %d | cells %d | reps per cell min %d max %d | lavaan %s\n\n", length(fs), length(reps),
            min(reps), max(reps), paste(unique(D$lavaan), collapse = ",")))

wilson <- function(x) { n <- length(x); if (n == 0) return(c(NA, NA, NA)); p <- mean(x); z <- 1.96
  c <- (p + z^2/(2*n)) / (1 + z^2/n); h <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n)
  100 * c(p, c - h, c + h) }
f <- function(v) !is.na(v) & v
clauses <- function(X, version) {
  w <- X$arm == "wlsmv_ord"
  cfi <- ifelse(w, if (version == "robust") X$cfi.rob else X$cfi.sc, X$cfi)
  tli <- ifelse(w, if (version == "robust") X$tli.rob else X$tli.sc, X$tli)
  rm  <- ifelse(w, if (version == "robust") X$rmsea.rob else X$rmsea.sc, X$rmsea)
  c_ <- f(cfi > .95); t_ <- f(tli > .95); r_ <- f(rm < .06); s_ <- f(X$srmr < .08)
  data.frame(c_, t_, r_, s_, cosmin = c_ | t_ | r_ | s_,
             rescue = r_ & !c_ & !t_ & !s_, srmr_sole = s_ & !c_ & !t_ & !r_,
             srmr_gt = f(X$srmr > rm))
}

# ---- fit quality ----
cat("=== Fit quality (percent of replications), by format x N, both arms ===\n")
q <- aggregate(cbind(ok = D$status == "ok", admiss = D$status == "ok" & D$admiss %in% TRUE,
                     degenerate = D$status == "degenerate", nonconv = D$status == "nonconv",
                     any_empty = D$n_empty > 0) ~ arm + fmt + n, data = D, FUN = function(x) round(100 * mean(x), 1))
print(q[order(q$arm, q$fmt, q$n), ], row.names = FALSE)

report <- function(E, version, label) {
  C <- clauses(E, version)
  cat(sprintf("\n==================== %s ====================\n", label))
  for (a in c("wlsmv_ord", "ml_ord_as_cont")) {
    g0 <- E$arm == a & E$kind == "none"
    cat(sprintf("\n--- %s: CORRECT models. Percent (95%% CI) ---\n", a))
    cat(sprintf("%-8s %2s %4s %4s | %-19s %-19s %-19s %-19s | %-19s | %-19s | %s\n", "format","p","load","N",
        "CFI fails","RMSEA fails","SRMR fails","COSMIN fails","RMSEA rescue","SRMR sole","SRMR>RMSEA"))
    for (ct in levels(E$fmt)) for (p in sort(unique(E$p))) for (ld in sort(unique(E$load))) for (n in sort(unique(E$n))) {
      g <- g0 & E$fmt == ct & E$p == p & E$load == ld & E$n == n; if (!any(g)) next
      s <- function(x) { w <- wilson(x); sprintf("%5.1f (%4.1f-%5.1f)", w[1], w[2], w[3]) }
      cat(sprintf("%-8s %2d %4.2f %4d | %s %s %s %s | %s | %s | %5.1f\n", ct, p, ld, n,
          s(!C$c_[g]), s(!C$r_[g]), s(!C$s_[g]), s(!C$cosmin[g]), s(C$rescue[g]), s(C$srmr_sole[g]),
          100 * mean(C$srmr_gt[g])))
    }
    g1 <- E$arm == a & E$conseq
    cat(sprintf("\n--- %s: CONSEQUENTIAL misspecification, COSMIN pass (false pass) %%, rows format x p x load x kind, cols N ---\n", a))
    tb <- tapply(C$cosmin[g1], list(paste(E$fmt[g1], sprintf("p%02d", E$p[g1]), E$load[g1], E$kind[g1]), E$n[g1]), mean)
    print(round(100 * tb, 1))
    gm <- E$arm == a & E$kind != "none" & !E$conseq
    if (any(gm)) { cat(sprintf("\n--- %s: non-consequential misfit, COSMIN pass %% ---\n", a))
      print(round(100 * tapply(C$cosmin[gm], list(paste(E$fmt[gm], sprintf("p%02d", E$p[gm]), E$load[gm], E$kind[gm]), E$n[gm]), mean), 1)) }
    cat(sprintf("\n--- %s: RMSEA rescue among ALL clean fits, and what it rescues (counts), by format x N ---\n", a))
    ga <- E$arm == a
    for (ct in levels(E$fmt)) for (n in sort(unique(E$n))) { g <- ga & E$fmt == ct & E$n == n
      cat(sprintf("  %-8s N=%3d  rescues %5d of %6d | correct %4d | consequential %4d | other misfit %4d\n", ct, n,
          sum(C$rescue[g]), sum(g), sum(C$rescue[g] & E$kind[g] == "none"), sum(C$rescue[g] & E$conseq[g]),
          sum(C$rescue[g] & E$kind[g] != "none" & !E$conseq[g]))) }
  }
  invisible(C)
}
E  <- D[D$status == "ok" & D$admiss %in% TRUE, ]
C1 <- report(E, "scaled", "PRIMARY: converged + admissible, WLSMV scaled indices")

# ---- pre-stated expectations ----
cat("\n=== Expectations (primary set, WLSMV scaled, correct models) ===\n")
g0 <- E$arm == "wlsmv_ord" & E$kind == "none"
cell <- interaction(E$fmt, E$p, E$load, E$n, drop = TRUE)
srmr_ff <- tapply(!C1$s_[g0], cell[g0], mean); resc <- tapply(C1$rescue[g0], cell[g0], mean)
parts <- strsplit(names(srmr_ff), ".", fixed = TRUE); fmt <- sapply(parts, `[`, 1); nn <- as.integer(sapply(parts, function(z) tail(z, 1)))
e1 <- fmt %in% c("c2_sym", "c2_skew", "c3_skew") & nn <= 200
cat(sprintf("E1 SRMR false-fail > 20%%: holds in %d of %d cells (binary sym/skew or 3-cat skew, N<=200); range %.1f-%.1f%%\n",
    sum(srmr_ff[e1] > .20), sum(e1), 100 * min(srmr_ff[e1]), 100 * max(srmr_ff[e1])))
cat(sprintf("E2 RMSEA rescue > 5%%: holds in %d of %d of those cells; range %.1f-%.1f%%\n",
    sum(resc[e1] > .05), sum(e1), 100 * min(resc[e1]), 100 * max(resc[e1])))
e3 <- fmt == "c5_sym"
cat(sprintf("E3 RMSEA rescue < 1%% for 5-cat symmetric: holds in %d of %d cells; max %.1f%%\n",
    sum(resc[e3] < .01), sum(e3), 100 * max(resc[e3])))

# ---- sensitivities ----
report(E[E$arm == "wlsmv_ord", ], "robust", "SENSITIVITY 1: WLSMV robust indices")
Eall <- D[D$status == "ok", ]
C3 <- clauses(Eall, "scaled"); g <- Eall$arm == "wlsmv_ord" & Eall$kind == "none"
cat("\n=== SENSITIVITY 2: all converged fits (incl. inadmissible), WLSMV scaled, correct models, by format x N ===\n")
print(round(100 * cbind(COSMIN_fail = tapply(!C3$cosmin[g], list(Eall$fmt[g], Eall$n[g]), mean),
                        rescue = tapply(C3$rescue[g], list(Eall$fmt[g], Eall$n[g]), mean)), 1))
saveRDS(list(D = D, E = E, C1 = C1), "out/study2b_analysis.rds")
