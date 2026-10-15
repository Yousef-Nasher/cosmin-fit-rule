# design2b.R -- Study 2b (item format) design: thresholds, conditions, population
# matrices, model syntax, job table. One-factor models are fitted throughout.
# 6 item formats x 3 scale lengths x 2 loadings x 3 model conditions (the two-factor
# condition only for 10 and 20 items) = 96 conditions; x 3 sample sizes = 288 cells.
CUTS <- list(c2_sym = c(0),              c2_skew = c(1.2),
             c3_sym = c(-0.75, 0.75),    c3_skew = c(0.4, 1.5),
             c5_sym = c(-1.5, -0.5, 0.5, 1.5), c5_skew = c(-0.2, 0.6, 1.2, 1.9))
PS    <- c(5, 10, 20)
LOADS <- c(0.45, 0.70)
NS    <- c(100, 200, 500)
MIS   <- list(list("none", 0, 0), list("resid", 1, .30), list("dim", 0, .70))
CHUNK <- 50                                  # replications per job
NCHUNK <- 4                                  # 200 replications per cell
SEED0 <- 20261001

build_pop <- function(p, load, kind, k, delta) {
  Lam <- matrix(load, p, 1); Phi <- matrix(1, 1, 1)
  if (kind == "dim") {                       # last half of items on a 2nd factor, r = delta
    kk <- p %/% 2; idx <- (p - kk + 1):p
    Lam <- cbind(Lam, 0); Lam[idx, 1] <- 0; Lam[idx, 2] <- load
    Phi <- matrix(c(1, delta, delta, 1), 2, 2)
  }
  common <- Lam %*% Phi %*% t(Lam); Sig <- common
  if (kind == "resid") { Sig[1, 2] <- Sig[2, 1] <- Sig[1, 2] + delta }
  diag(Sig) <- 1
  nm <- paste0("x", 1:p); dimnames(Sig) <- list(nm, nm)
  list(Sig = Sig, common = common)
}
model_syntax <- function(p) paste0("F1 =~ ", paste0("x", 1:p, collapse = " + "))

conditions <- function() {
  out <- list()
  for (ct in names(CUTS)) for (p in PS) for (load in LOADS) for (ms in MIS) {
    kind <- ms[[1]]; k <- ms[[2]]; delta <- ms[[3]]
    if (kind == "dim" && p < 10) next
    id <- sprintf("%s_p%02d_l%.2f_%s", ct, p, load, kind)
    out[[id]] <- list(id = id, ct = ct, p = p, load = load, kind = kind, k = k, delta = delta)
  }
  out
}
job_table <- function(n_chunks = NCHUNK) {
  C <- conditions(); cells <- list(); i <- 0
  for (c in C) for (n in NS) { i <- i + 1
    cells[[i]] <- c(c, list(n = n, cell_index = i, cell = sprintf("%s_n%04d", c$id, n))) }
  J <- list(); j <- 0
  for (ch in 1:n_chunks) for (cl in cells) { j <- j + 1
    J[[j]] <- c(cl, list(chunk = ch, job = sprintf("%s_c%02d", cl$cell, ch),
                         seed = SEED0 + 1000 * ch + cl$cell_index)) }
  J
}
