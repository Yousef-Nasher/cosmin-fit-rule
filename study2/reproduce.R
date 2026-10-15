# reproduce.R -- Study 2: re-run the analysis from the frozen simulation output and compare
# every result with the stored copy.
# Run from study2/:  Rscript reproduce.R        (base R only; about one minute)
# The independent Python check (pandas, numpy) is skipped with a note if Python is not available.
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, out = NULL) {
  status <- system2(rscript, script, stdout = if (is.null(out)) "" else out)
  if (status != 0) stop(script, " failed", call. = FALSE)
}
unlink("out", recursive = TRUE); dir.create("out")
unzip("data/study2_final_chunks.zip", exdir = "out")

run("code/analysis.R",               "out/study2_results.txt")
run("code/analysis_extra.R",         "out/study2_extra.txt")
run("code/hu_bentler_sensitivity.R", "out/hu_bentler_sensitivity.txt")
run("code/table_rules.R",            "out/table_rules.log")

same_text <- function(f) identical(readLines(file.path("out", f)), readLines(file.path("results", f)))
same_csv  <- function(f) isTRUE(all.equal(read.csv(file.path("out", f)), read.csv(file.path("results", f)),
                                          tolerance = 1e-10))
ok <- c(study2_results.txt = same_text("study2_results.txt"),
        study2_extra.txt = same_text("study2_extra.txt"),
        hu_bentler_sensitivity.txt = same_text("hu_bentler_sensitivity.txt"),
        table_rules_by_N.csv = same_csv("table_rules_by_N.csv"))
py <- Filter(function(p) nzchar(p) && system2(p, c("-c", shQuote("import pandas, numpy")),
                                                stdout = FALSE, stderr = FALSE) == 0,
             Sys.which(c("python3", "python")))
if (length(py)) {
  if (system2(py[[1]], "code/verify_independent.py", stdout = "out/verify_independent.txt") != 0)
    stop("code/verify_independent.py failed", call. = FALSE)
  ok["verify_independent.txt"] <- same_text("verify_independent.txt")
} else cat("Python with pandas and numpy not found: independent check skipped\n")
for (f in names(ok)) cat(sprintf("%-28s %s\n", f, if (ok[[f]]) "identical to stored copy" else "DIFFERS from stored copy"))
if (!all(ok)) quit(status = 1)
