# worker2b.R -- Study 2b simulation worker.
# Run from study2b/code/:  Rscript worker2b.R <worker_name> [rev]
# Jobs are claimed by lock folders and written one CSV per job to out/chunks/, as in Study 2.
source("design2b.R"); source("simulate2b.R")
a <- commandArgs(TRUE); WNAME <- a[1]
OUT <- "out/chunks"; LCK <- "out/locks"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(LCK, recursive = TRUE, showWarnings = FALSE)

JOBS <- job_table()
if (length(a) >= 2 && a[2] == "rev") JOBS <- rev(JOBS)
cat(sprintf("[%s] %d jobs in table, lavaan %s, started %s\n", WNAME, length(JOBS),
    packageVersion("lavaan"), format(Sys.time(), "%H:%M:%S")))
for (J in JOBS) {
  fout <- file.path(OUT, paste0(J$job, ".csv"))
  if (file.exists(fout)) next
  if (!dir.create(file.path(LCK, J$job), showWarnings = FALSE)) next
  t0 <- Sys.time(); D <- run_job(J)
  D$seed <- J$seed; D$worker <- WNAME; D$lavaan <- as.character(packageVersion("lavaan"))
  tmp <- paste0(fout, ".tmp"); write.csv(D, tmp, row.names = FALSE); file.rename(tmp, fout)
  cat(sprintf("[%s] %s  %4.0fs  %s\n", WNAME, J$job,
      as.numeric(difftime(Sys.time(), t0, units = "secs")), format(Sys.time(), "%H:%M")))
}
cat(sprintf("[%s] no jobs left. finished %s\n", WNAME, format(Sys.time(), "%H:%M:%S")))
