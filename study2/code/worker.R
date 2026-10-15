# worker.R -- Study 2 simulation worker.
# Run from study2/code/:  Rscript worker.R <worker_name> [n_chunks] [rev]
# Each job (one cell x 50 replications) is claimed by creating a lock folder and written
# as one CSV to out/chunks/, through a temporary file. Any number of workers can therefore
# run, stop and restart side by side; "rev" makes a worker take the jobs in reverse order.
source("design.R"); source("simulate.R")
a <- commandArgs(TRUE); WNAME <- a[1]
NCH <- if (length(a) >= 2) as.integer(a[2]) else NCHUNK
OUT <- "out/chunks"; LCK <- "out/locks"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(LCK, recursive = TRUE, showWarnings = FALSE)

JOBS <- job_table(NCH)
if (length(a) >= 3 && a[3] == "rev") JOBS <- rev(JOBS)
cat(sprintf("[%s] %d jobs in table, lavaan %s, started %s\n", WNAME, length(JOBS),
    packageVersion("lavaan"), format(Sys.time(), "%H:%M:%S")))
for (J in JOBS) {
  fout <- file.path(OUT, paste0(J$job, ".csv"))
  if (file.exists(fout)) next
  if (!dir.create(file.path(LCK, J$job), showWarnings = FALSE)) next   # claimed by another
  t0 <- Sys.time()
  D <- run_job(J)
  D$seed <- J$seed; D$worker <- WNAME
  tmp <- paste0(fout, ".tmp"); write.csv(D, tmp, row.names = FALSE); file.rename(tmp, fout)
  cat(sprintf("[%s] %s  %4.0fs  %s\n", WNAME, J$job,
      as.numeric(difftime(Sys.time(), t0, units = "secs")), format(Sys.time(), "%H:%M")))
}
cat(sprintf("[%s] no jobs left. finished %s\n", WNAME, format(Sys.time(), "%H:%M:%S")))
