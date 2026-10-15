# reproduce.R -- Study 3: re-run the draw and the analysis from the stored data and compare every
# result with the stored copy.
# Run from study3/:  Rscript reproduce.R        (needs the cosminfit package; a few seconds)
# The draw check and the independent check need Python 3 (the latter with pandas and numpy);
# they are skipped with a note if Python is not available.
rscript <- file.path(R.home("bin"), "Rscript")
unlink("out", recursive = TRUE); dir.create("out")
if (system2(rscript, c("code/analysis.R", "out"), stdout = FALSE) != 0) stop("code/analysis.R failed", call. = FALSE)
if (system2(rscript, c("code/exploratory_q1_period.R", "out"), stdout = FALSE) != 0)
  stop("code/exploratory_q1_period.R failed", call. = FALSE)

same_text <- function(f) identical(readLines(file.path("out", f)), readLines(file.path("results", f)))
same_csv <- function(f) isTRUE(all.equal(read.csv(file.path("out", f)), read.csv(file.path("results", f)),
                                         tolerance = 1e-10))
stored <- list.files("results")
ok <- sapply(stored, function(f) {
  if (!file.exists(file.path("out", f))) return(FALSE)
  if (grepl("[.]csv$", f)) same_csv(f) else if (f == "verify_independent.txt") TRUE else same_text(f)
})

python <- function(modules) {
  py <- Filter(function(p) nzchar(p) && system2(p, c("-c", shQuote(paste("import", modules))),
                                                  stdout = FALSE, stderr = FALSE) == 0,
               Sys.which(c("python3", "python")))
  if (length(py)) py[[1]] else ""
}
py <- python("csv")
if (nzchar(py)) {
  if (system2(py, c("code/draw_sample.py", "out/draw.csv")) != 0) stop("code/draw_sample.py failed", call. = FALSE)
  d <- read.csv("out/draw.csv", stringsAsFactors = FALSE)
  s <- read.csv("data/study3_fulltext_screening.csv", stringsAsFactors = FALSE)
  m <- merge(s, d, by = "frame_id", suffixes = c("", ".draw"))
  ok["draw order (code/draw_sample.py)"] <- nrow(m) == nrow(s) && nrow(d) == nrow(s) &&
    all(m$year_stratum == m$year_stratum.draw) && all(m$draw_list == m$draw_list.draw) && all(m$rank == m$rank.draw)
} else cat("Python not found: draw check skipped\n")
py <- python("pandas, numpy")
if (nzchar(py)) {
  if (system2(py, "code/verify_independent.py", stdout = "out/verify_independent.txt") != 0)
    stop("code/verify_independent.py failed", call. = FALSE)
  ok["verify_independent.txt"] <- same_text("verify_independent.txt")
} else {
  ok <- ok[names(ok) != "verify_independent.txt"]
  cat("Python with pandas and numpy not found: independent check skipped\n")
}
for (f in names(ok)) cat(sprintf("%-44s %s\n", f, if (ok[[f]]) "identical to stored copy" else "DIFFERS from stored copy"))
if (!all(ok)) quit(status = 1)
