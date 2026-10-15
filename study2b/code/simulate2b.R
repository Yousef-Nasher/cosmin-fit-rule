# simulate2b.R -- Study 2b: generate the data of one job and fit the two analysis arms.
# Sourced by worker2b.R and rerun_job2b.R (after design2b.R).
# Latent responses are cut at the thresholds of the item format. Each data set is analysed
# with WLSMV on the ordinal items and with ML on the items treated as continuous.
# Categories that are not observed are left as they fall; a data set with a constant item
# is marked "degenerate" and not fitted.
suppressPackageStartupMessages(library(lavaan))

EMPTY <- list(status=NA_character_, admiss=NA, chisq=NA_real_, df=NA_real_,
  cfi=NA_real_, tli=NA_real_, rmsea=NA_real_, srmr=NA_real_,
  chisq.sc=NA_real_, cfi.sc=NA_real_, tli.sc=NA_real_, rmsea.sc=NA_real_,
  cfi.rob=NA_real_, tli.rob=NA_real_, rmsea.rob=NA_real_)
harvest <- function(f) {
  o <- EMPTY
  if (inherits(f, "error"))                 { o$status <- "error";   return(o) }
  if (!isTRUE(lavInspect(f, "converged")))  { o$status <- "nonconv"; return(o) }
  sd <- tryCatch(lavInspect(f, "std"), error = function(e) NULL)
  if (!is.null(sd)) o$admiss <- all(abs(sd$lambda) < 0.999) && all(diag(sd$theta) > -1e-8)
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
  B <- build_pop(J$p, J$load, J$kind, J$k, J$delta)
  L <- chol(B$Sig); syn <- model_syntax(J$p); p <- J$p; n <- J$n; cuts <- CUTS[[J$ct]]
  nm <- paste0("x", 1:p)
  meta <- function(rep, arm, ne) list(cell=J$cell, ct=J$ct, p=p, load=J$load, kind=J$kind,
    k=J$k, delta=J$delta, n=n, chunk=J$chunk, rep=rep, arm=arm, n_empty=ne)
  rows <- list(); r <- 0
  for (rep in (J$chunk - 1) * CHUNK + 1:CHUNK) {
    Z  <- matrix(rnorm(n * p), n, p) %*% L
    Xi <- apply(Z, 2, function(v) findInterval(v, cuts)); colnames(Xi) <- nm
    nobs <- apply(Xi, 2, function(v) length(unique(v)))
    ne   <- sum(nobs < length(cuts) + 1)          # items with an unobserved category
    degenerate <- any(nobs < 2)
    Xo <- as.data.frame(lapply(as.data.frame(Xi), ordered))
    h <- if (degenerate) modifyList(EMPTY, list(status = "degenerate")) else
         harvest(fit(cfa(syn, data = Xo, ordered = nm, estimator = "WLSMV", std.lv = TRUE)))
    r <- r + 1; rows[[r]] <- c(meta(rep, "wlsmv_ord", ne), h)
    h <- if (degenerate) modifyList(EMPTY, list(status = "degenerate")) else
         harvest(fit(cfa(syn, data = as.data.frame(Xi), std.lv = TRUE)))
    r <- r + 1; rows[[r]] <- c(meta(rep, "ml_ord_as_cont", ne), h)
  }
  do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
}
