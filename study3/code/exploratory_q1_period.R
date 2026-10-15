# Study 3, exploratory analysis that was not registered: share of reviews that state the structural
# validity criterion, by publication period (2018 to 2024, before PRISMA-COSMIN for OMIs 2024 could
# affect most reviews, and 2025 to 2026).
# Run from study3/:  Rscript code/exploratory_q1_period.R out
# Reads data/study3_reviews.csv and writes exploratory_q1_period.txt to the folder given as the first
# argument (default "results").
suppressPackageStartupMessages(library(cosminfit))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "results"
dir.create(out_dir, showWarnings = FALSE)

reviews <- read.csv("data/study3_reviews.csv", stringsAsFactors = FALSE, na.strings = "")
reviews$stated <- reviews$wording != "not stated"
reviews$period <- ifelse(reviews$year <= 2024, "2018-2024", "2025-2026")

pct <- function(x) sprintf("%.1f%%", 100 * x)
out <- c("Exploratory (not registered): reviews that state the criterion, by publication period",
         "Stated = wording coded disjunctive, conjunctive or mixed; Wilson 95% CI", "")
for (p in c("2018-2024", "2025-2026")) {
  s <- reviews[reviews$period == p, ]
  k <- sum(s$stated)
  w <- wilson_ci(k, nrow(s))
  out <- c(out, sprintf("  %s  %2d/%2d = %s (95%% CI %s to %s)", p, k, nrow(s), pct(w$est),
                        pct(w$lower), pct(w$upper)))
}
tab <- table(reviews$period, reviews$stated)
out <- c(out, "", sprintf("  Fisher's exact test: p = %.2f", fisher.test(tab)$p.value))
writeLines(out, file.path(out_dir, "exploratory_q1_period.txt"))
