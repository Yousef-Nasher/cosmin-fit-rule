# design.R -- Study 2 design: conditions, population matrices, model syntax, job table.
# 4 scale shapes x 2 loadings x 6 model conditions, minus the combinations that do not
# apply (see conditions()), gives 38 conditions; x 4 sample sizes = 152 cells.
SHAPES <- list(c(5,1), c(10,1), c(10,2), c(12,2))   # (items, factors)
LOADS  <- c(0.50, 0.70)
NS     <- c(100, 200, 500, 1000)
# (kind, k, delta): correct model; k residual correlations of delta; k cross-loadings of
# delta; an omitted dimension correlated delta with its parent factor
MIS    <- list(list("none",0,0), list("resid",1,.30), list("resid",2,.30),
               list("resid",3,.40), list("cross",2,.40), list("dim",0,.70))
ORD_CUTS <- c(-1.5, -0.5, 0.5, 1.5)   # 5 symmetric categories
CHUNK  <- 50                          # replications per job
NCHUNK <- 10                          # 500 replications per cell

simple_pattern <- function(p, m) {
  pat <- matrix(0, p, m); per <- p %/% m
  for (j in 1:m) { lo <- (j-1)*per + 1; hi <- if (j < m) j*per else p
                   pat[lo:hi, j] <- 1 }
  pat
}

# returns population correlation matrix AND the true generating parameters
build_pop <- function(p, m, load, phi, kind, k, delta) {
  pat <- simple_pattern(p, m); Lam <- pat * load
  Phi <- matrix(phi, m, m); diag(Phi) <- 1
  if (kind == "cross") { src <- which(pat[, 2] == 1)[1:k]; Lam[src, 1] <- delta }
  if (kind == "dim") {
    kk  <- max(2, (p %/% m) %/% 2); idx <- tail(which(pat[, 1] == 1), kk)
    Lam <- cbind(Lam, 0); Lam[idx, 1] <- 0; Lam[idx, m+1] <- load
    Phi2 <- matrix(phi, m+1, m+1); diag(Phi2) <- 1
    Phi2[1, m+1] <- Phi2[m+1, 1] <- delta; Phi <- Phi2
  }
  common <- Lam %*% Phi %*% t(Lam)
  Sig <- common
  if (kind == "resid") {
    prs <- list(c(1,2), c(3,4), c(5,6))[1:k]
    for (pr in prs) { Sig[pr[1], pr[2]] <- Sig[pr[1], pr[2]] + delta
                      Sig[pr[2], pr[1]] <- Sig[pr[1], pr[2]] }
  }
  diag(Sig) <- 1
  nm <- paste0("x", 1:p); dimnames(Sig) <- list(nm, nm)
  list(Sig = Sig, common = common, Lam = Lam, Phi = Phi, pat = pat)
}

model_syntax <- function(p, m) {
  pat <- simple_pattern(p, m)
  paste(sapply(1:m, function(j)
    paste0("F", j, " =~ ", paste0("x", which(pat[, j] == 1), collapse = " + "))),
    collapse = "\n")
}

conditions <- function() {
  out <- list()
  for (sh in SHAPES) for (load in LOADS) for (ms in MIS) {
    p <- sh[1]; m <- sh[2]; kind <- ms[[1]]; k <- ms[[2]]; delta <- ms[[3]]
    if (kind == "cross" && m < 2) next
    if (kind == "dim"   && p %/% m < 6) next
    if (kind == "resid" && 2*k > p) next
    phi <- if (m > 1) 0.30 else 0
    B <- build_pop(p, m, load, phi, kind, k, delta)
    if (min(eigen(B$Sig, only.values = TRUE)$values) <= 1e-6) next
    id <- sprintf("p%02d_m%d_l%.2f_%s_k%d_d%.2f", p, m, load, kind, k, delta)
    out[[id]] <- list(id=id, p=p, m=m, load=load, phi=phi, kind=kind, k=k, delta=delta)
  }
  out
}

# Job list and seeds. Cell order is fixed by conditions() x NS. One job = one cell x one
# chunk of 50 replications. seed = 20260927 + 1000*chunk + cell_index, unique per job, so
# any job can be re-run on its own and the seed of a job does not depend on n_chunks.
job_table <- function(n_chunks = NCHUNK) {
  C <- conditions(); cells <- list(); i <- 0
  for (c in C) for (n in NS) { i <- i + 1
    cells[[i]] <- c(c, list(n = n, cell_index = i,
                            cell = sprintf("%s_n%04d", c$id, n))) }
  J <- list(); j <- 0
  for (ch in 1:n_chunks) for (cl in cells) { j <- j + 1
    J[[j]] <- c(cl, list(chunk = ch, job = sprintf("%s_c%02d", cl$cell, ch),
                         seed = 20260927 + 1000*ch + cl$cell_index)) }
  J
}
