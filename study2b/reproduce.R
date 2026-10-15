# reproduce.R -- Study 2b: re-run the analysis from the frozen simulation output and compare
# the result with the stored copy.
# Run from study2b/:  Rscript reproduce.R        (base R only; about one minute)
# The independent Python check (pandas) is skipped with a note if Python is not available.
rscript <- file.path(R.home("bin"), "Rscript")
unlink("out", recursive = TRUE); dir.create("out")
unzip("data/study2b_chunks.zip", exdir = "out")
status <- system2(rscript, "code/analysis2b.R", stdout = "out/study2b_results.txt")
if (status != 0) stop("code/analysis2b.R failed", call. = FALSE)
same_text <- function(f) identical(readLines(file.path("out", f)), readLines(file.path("results", f)))
ok <- c(study2b_results.txt = same_text("study2b_results.txt"))
py <- Filter(function(p) nzchar(p) && system2(p, c("-c", shQuote("import pandas")),
                                                stdout = FALSE, stderr = FALSE) == 0,
             Sys.which(c("python3", "python")))
if (length(py)) {
  if (system2(py[[1]], "code/verify_independent_2b.py", stdout = "out/verify_independent_2b.txt") != 0)
    stop("code/verify_independent_2b.py failed", call. = FALSE)
  ok["verify_independent_2b.txt"] <- same_text("verify_independent_2b.txt")
} else cat("Python with pandas not found: independent check skipped\n")
for (f in names(ok)) cat(sprintf("%-28s %s\n", f, if (ok[[f]]) "identical to stored copy" else "DIFFERS from stored copy"))
if (!all(ok)) quit(status = 1)
