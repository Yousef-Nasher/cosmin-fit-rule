# rerun_job2b.R -- Study 2b: re-run one simulation job from its seed and compare it with
# the frozen output.
# Run from study2b/code/:  Rscript rerun_job2b.R <job>
#   e.g. Rscript rerun_job2b.R c2_skew_p05_l0.45_none_n0100_c01
source("design2b.R"); source("simulate2b.R")
job <- commandArgs(TRUE)[1]
J <- Filter(function(j) j$job == job, job_table())
if (length(J) != 1) stop("unknown job: ", job)
new <- run_job(J[[1]])
old <- read.csv(unz("../data/study2b_chunks.zip", paste0("chunks/", job, ".csv")), stringsAsFactors = FALSE)
num <- c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "chisq.sc", "cfi.sc", "tli.sc", "rmsea.sc",
         "cfi.rob", "tli.rob", "rmsea.rob")
same_status <- identical(as.character(new$status), as.character(old$status))
both <- new$status == "ok" & old$status == "ok"
maxdiff <- max(sapply(num, function(v) suppressWarnings(max(abs(new[[v]][both] - old[[v]][both]), na.rm = TRUE))))
cut <- function(d) cbind(d$cfi > .95, d$tli > .95, d$rmsea < .06, d$srmr < .08, d$cfi.sc > .95, d$tli.sc > .95,
                         d$rmsea.sc < .06, d$cfi.rob > .95, d$rmsea.rob < .06)[both, , drop = FALSE]
flips <- sum(cut(new) != cut(old), na.rm = TRUE)
cat(sprintf("job %s | lavaan %s here, %s in the frozen file\n", job, packageVersion("lavaan"), old$lavaan[1]))
cat(sprintf("rows %d | fit status identical: %s | largest difference in any index: %.3g | cut-off decisions that differ: %d\n",
            nrow(new), same_status, maxdiff, flips))
if (!same_status || flips > 0 || maxdiff > 1e-4) quit(status = 1)
