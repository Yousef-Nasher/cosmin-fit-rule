# Base-R test runner: the tests need no package beyond cosminfit itself.
# Each file in tests/unit/ uses check(), check_error() and check_warning().
library(cosminfit)

n_passed <- 0
check <- function(cond, label) {
  if (!isTRUE(all(cond))) stop("Test failed: ", label, call. = FALSE)
  n_passed <<- n_passed + 1
  invisible(TRUE)
}
check_error <- function(expr, label) {
  failed <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  check(failed, paste("expected error:", label))
}
check_warning <- function(expr, label) {
  warned <- FALSE
  withCallingHandlers(force(expr), warning = function(w) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  check(warned, paste("expected warning:", label))
}
check_no_warning <- function(expr, label) {
  warned <- FALSE
  withCallingHandlers(force(expr), warning = function(w) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  check(!warned, paste("unexpected warning:", label))
}

# the unit tests sit next to this script; find them from wherever the script is run
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]))) else getwd()
files <- sort(list.files(file.path(here, "unit"), pattern = "^test-.*[.]R$", full.names = TRUE))
if (!length(files)) stop("No test files found in ", file.path(here, "unit"), call. = FALSE)
for (f in files) {
  before <- n_passed
  sys.source(f, envir = environment())
  cat(sprintf("%-32s %4d checks passed\n", basename(f), n_passed - before))
}
cat(sprintf("TOTAL: %d checks passed, 0 failed\n", n_passed))
