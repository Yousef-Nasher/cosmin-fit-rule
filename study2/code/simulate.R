# simulate.R -- Study 2: generate the data of one job and fit the three analysis arms.
# Sourced by worker.R and rerun_job.R (after design.R).
# Every data set is analysed three ways: ML on the continuous items, ML on the
# five-category items treated as continuous, and WLSMV on the five-category items.
suppressPackageStartupMessages(library(lavaan))

EMPTY <- list(status=NA_character_, admiss=NA, chisq=NA_real_, df=NA_real_,
  cfi=NA_real_, tli=NA_real_, rmsea=NA_real_, srmr=NA_real_,
  chisq.sc=NA_real_, cfi.sc=NA_real_, tli.sc=NA_real_, rmsea.sc=NA_real_,
  cfi.rob=NA_real_, tli.rob=NA_real_, rmsea.rob=NA_real_)

harvest <- function(f) {
  o <- EMPTY
  if (inherits(f, "error"))            { o$status <- "error";   return(o) }
  if (!isTRUE(lavInspect(f, "converged"))) { o$status <- "nonconv"; return(o) }
  sd <- tryCatch(lavInspect(f, "std"), error = function(e) NULL)
  if (!is.null(sd)) {
    ok <- all(abs(sd$lambda) < 0.999) && all(diag(sd$theta) > -1e-8)
    if (ncol(sd$psi) > 1) ok <- ok && all(abs(sd$psi[lower.tri(sd$psi)]) < 0.999)
    o$admiss <- ok
  }
  fm <- tryCatch(fitMeasures(f), error = function(e) NULL)
  if (is.null(fm)) { o$status <- "fm_error"; return(o) }
  g <- function(k) if (k %in% names(fm)) as.numeric(fm[[k]]) else NA_real_
  o$chisq <- g("chisq"); o$df <- g("df"); o$cfi <- g("cfi"); o$tli <- g("tli")
  o$rmsea <- g("rmsea"); o$srmr <- g("srmr")
  o$chisq.sc <- g("chisq.scaled"); o$cfi.sc <- g("cfi.scaled")
  o$tli.sc <- g("tli.scaled"); o$rmsea.sc <- g("rmsea.scaled")
  o$cfi.rob <- g("cfi.robust"); o$tli.rob <- g("tli.robust"); o$rmsea.rob <- g("rmsea.robust")
  o$status <- "ok"; o
}
fit <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)

run_job <- function(J) {
  set.seed(J$seed)
  B <- build_pop(J$p, J$m, J$load, J$phi, J$kind, J$k, J$delta)
  L <- chol(B$Sig); syn <- model_syntax(J$p, J$m); p <- J$p; n <- J$n
  nm <- paste0("x", 1:p)
  meta <- function(rep, arm) list(cell=J$cell, p=p, m=J$m, load=J$load,
    kind=J$kind, k=J$k, delta=J$delta, n=n, chunk=J$chunk, rep=rep, arm=arm)
  rows <- list(); r <- 0
  for (rep in (J$chunk - 1) * CHUNK + 1:CHUNK) {
    Z <- matrix(rnorm(n * p), n, p) %*% L; colnames(Z) <- nm
    # arm 1: continuous indicators, ML on Pearson
    h <- harvest(fit(cfa(syn, data = as.data.frame(Z), std.lv = TRUE)))
    r <- r + 1; rows[[r]] <- c(meta(rep, "ml_cont"), h)
    # same respondents, 5-category ordinal items
    Xi <- apply(Z, 2, function(v) findInterval(v, ORD_CUTS)); colnames(Xi) <- nm
    degenerate <- any(apply(Xi, 2, function(v) length(unique(v))) < 2)
    # arm 2: ordinal items wrongly treated as continuous, ML
    h <- if (degenerate) modifyList(EMPTY, list(status = "degenerate")) else
         harvest(fit(cfa(syn, data = as.data.frame(Xi), std.lv = TRUE)))
    r <- r + 1; rows[[r]] <- c(meta(rep, "ml_ord_as_cont"), h)
    # arm 3: ordinal items, WLSMV on polychoric correlations
    if (degenerate) h <- modifyList(EMPTY, list(status = "degenerate")) else {
      Xo <- as.data.frame(lapply(as.data.frame(Xi), ordered))
      h <- harvest(fit(cfa(syn, data = Xo, ordered = nm, estimator = "WLSMV",
                           std.lv = TRUE)))
    }
    r <- r + 1; rows[[r]] <- c(meta(rep, "wlsmv_ord"), h)
  }
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}
