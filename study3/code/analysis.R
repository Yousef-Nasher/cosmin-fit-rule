# Study 3: how published COSMIN reviews apply the structural validity criterion.
# Run from study3/:  Rscript code/analysis.R out
# Reads data/, writes the results to the folder given as the first argument (default "results").
suppressPackageStartupMessages(library(cosminfit))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "results"
dir.create(out_dir, showWarnings = FALSE)
SEED <- 545749
B <- 2000
RULES <- c("cosmin", "no_rmsea", "srmr_df", "two_of_three")
REPAIRED <- RULES[-1]

reviews <- read.csv("data/study3_reviews.csv", stringsAsFactors = FALSE, na.strings = "")
cfa <- read.csv("data/study3_cfa_results.csv", stringsAsFactors = FALSE, na.strings = "")
inst <- read.csv("data/study3_instruments.csv", stringsAsFactors = FALSE, na.strings = "")
for (v in c("value_flags", "rating", "rating_scope", "instrument_group"))
  cfa[[v]][is.na(cfa[[v]])] <- ""
for (v in c("sv", "ic", "recommendation", "sv_quality"))
  inst[[v]][is.na(inst[[v]])] <- ""

report <- character()
say <- function(...) report <<- c(report, sprintf(...))
pct <- function(x) sprintf("%.1f%%", 100 * x)
ci_text <- function(est, lo, hi) sprintf("%s (95%% CI %s to %s)", pct(est), pct(lo), pct(hi))
boot <- function(x, cluster) {
  b <- cluster_boot_prop(x, cluster, B = B, seed = SEED)
  b$text <- sprintf("%d/%d = %s, %d reviews", sum(x, na.rm = TRUE), b$n,
                    ci_text(b$est, b$lower, b$upper), b$n_clusters)
  b
}
write_out <- function(d, f) write.csv(d, file.path(out_dir, f), row.names = FALSE, na = "")

has_flag <- function(flags, what) grepl(what, flags, fixed = TRUE)
cfa$range <- has_flag(cfa$value_flags, "range")
cfa$implausible <- has_flag(cfa$value_flags, "implausible")

# ---- Sample -------------------------------------------------------------------------------
screen <- read.csv("data/study3_fulltext_screening.csv", stringsAsFactors = FALSE, na.strings = "")
say("STUDY 3 RESULTS")
say("")
say("Sample")
say("  Frame: %d records; drawn: %d; replacement queue: %d",
    nrow(screen), sum(screen$draw_list == "drawn"), sum(screen$draw_list == "replacement queue"))
tab <- table(factor(screen$fulltext_status,
                    c("included", "excluded", "full text not obtained", "not reached")))
for (s in names(tab)) say("  %-24s %4d", s, tab[[s]])
ex <- sort(table(screen$reason[screen$fulltext_status == "excluded"]), decreasing = TRUE)
for (s in names(ex)) say("    excluded: %-82s %3d", s, ex[[s]])
say("  Included reviews: %d (drawn %d, replacement %d); COSMIN version 1: %d, version 2: %d",
    nrow(reviews), sum(reviews$entered_as == "drawn"), sum(reviews$entered_as == "replacement"),
    sum(reviews$cosmin_version == 1), sum(reviews$cosmin_version == 2))
say("  CFA rows: %d (study level %d, pooled %d, summary %d); instruments with a review-level record: %d",
    nrow(cfa), sum(cfa$level == "study"), sum(cfa$level == "pooled"), sum(cfa$level == "summary"),
    nrow(inst))
strata <- table(reviews$year_stratum)
say("  Included by year stratum: %s", paste(names(strata), strata, sep = " ", collapse = ", "))
search <- read.csv("data/study3_search_counts.csv", stringsAsFactors = FALSE)
reached <- screen$fulltext_status != "not reached"
flow <- rbind(
  data.frame(stage = search$stage, detail = search$source, n = search$n),
  data.frame(stage = "records reached in draw order", detail = "frame", n = sum(reached)),
  data.frame(stage = "records not reached (stratum quota filled)", detail = "frame", n = sum(!reached)),
  data.frame(stage = "full text not obtained", detail = "reached", n = tab[["full text not obtained"]]),
  data.frame(stage = "full text assessed", detail = "reached", n = tab[["excluded"]] + tab[["included"]]),
  data.frame(stage = paste("excluded:", names(ex)), detail = "full text", n = as.vector(ex)),
  data.frame(stage = "reviews included", detail = "full text", n = tab[["included"]]),
  data.frame(stage = "study-level CFA results", detail = "included reviews", n = sum(cfa$level == "study")))
stopifnot(sum(reached) == sum(tab[c("included", "excluded", "full text not obtained")]),
          search$n[search$stage == "records in the sampling frame"] == nrow(screen))
write_out(flow, "flow_counts.csv")

# ---- Q1 ---------------------------------------------------------------------------------
say("")
say("Q1. How reviews state the criterion (unit: review, n = %d; Wilson 95%% CI)", nrow(reviews))
lev <- c("disjunctive", "conjunctive", "mixed", "availability", "not stated")
k <- as.vector(table(factor(reviews$wording, lev)))
q1 <- data.frame(wording = lev, n = k, wilson_ci(k, nrow(reviews)))
for (i in seq_len(nrow(q1)))
  say("  %-13s %3d  %s", q1$wording[i], q1$n[i], ci_text(q1$est[i], q1$lower[i], q1$upper[i]))
unamb <- sum(reviews$wording %in% c("disjunctive", "conjunctive"))
w <- wilson_ci(unamb, nrow(reviews))
say("  Stated unambiguously (disjunctive or conjunctive): %d/%d = %s", unamb, nrow(reviews),
    ci_text(w$est, w$lower, w$upper))
q1_by_version <- table(version = reviews$cosmin_version, wording = factor(reviews$wording, lev))
for (v in rownames(q1_by_version))
  say("  COSMIN version %s: %s", v, paste(colnames(q1_by_version), q1_by_version[v, ], sep = " ", collapse = ", "))
write_out(q1, "q1_wording.csv")

# ---- Results used for Q2 to Q4 -------------------------------------------------------------
# One CFA result = a study-level row with single values (rows that print a range summarise
# several models and are kept for a sensitivity analysis only).
res <- cfa[cfa$level == "study" & !cfa$range, ]
say("")
say("CFA results analysed: %d study-level rows from %d reviews (%d rows with a printed range set aside)",
    nrow(res), length(unique(res$review_id)), sum(cfa$level == "study" & cfa$range))

# degrees of freedom for the df-scaled SRMR rule: as reported, otherwise the one-factor model on
# the reported number of items (the largest df any model on those items can have, so the most
# lenient scaled cut-off)
prep <- function(d) {
  d$p <- suppressWarnings(as.numeric(d$items))
  d$df_used <- d$df
  one <- is.na(d$df_used) & !is.na(d$p) & d$p >= 4
  d$df_used[one] <- d$p[one] * (d$p[one] - 3) / 2
  d
}
rate_all <- function(d, missing = "reported") {
  d <- prep(d)
  r <- rerate(d, rules = RULES, missing = missing,
              cols = list(cfi = "cfi", tli = "tli", rmsea = "rmsea", srmr = "srmr", p = "p", df = "df_used"))
  # the scaled SRMR term is unknown, not absent, when SRMR is reported without p or df:
  # an insufficient rating that rests on this gap is indeterminate
  gap <- !is.na(r$srmr) & r$srmr < 0.08 & (is.na(r$p) | is.na(r$df_used))
  fix <- gap & r$rating_srmr_df == "insufficient"
  r$rating_srmr_df[fix] <- "indeterminate"
  r
}
res <- rate_all(res)

# ---- Q2 ---------------------------------------------------------------------------------
say("")
say("Q2. Number of the four COSMIN indices (CFI, TLI, RMSEA, SRMR) reported per CFA result (n = %d)", nrow(res))
q2 <- as.data.frame(table(n_indices = factor(res$n_indices, 0:4)), responseName = "n")
q2$percent <- round(100 * q2$n / nrow(res), 1)
for (i in seq_len(nrow(q2))) say("  %s indices: %4d  (%5.1f%%)", q2$n_indices[i], q2$n[i], q2$percent[i])
qs <- quantile(res$n_indices, c(.25, .5, .75), type = 2)
say("  Median %g (IQR %g to %g); mean %.2f", qs[2], qs[1], qs[3], mean(res$n_indices))
say("  Reported: CFI %d, TLI %d, RMSEA %d, SRMR %d",
    sum(!is.na(res$cfi)), sum(!is.na(res$tli)), sum(!is.na(res$rmsea)), sum(!is.na(res$srmr)))
two <- boot(res$n_indices >= 2, res$review_id)
say("  Two or more indices: %s", two$text)
say("  All three of CFI, RMSEA and SRMR: %d; items reported: %d; df reported: %d",
    sum(!is.na(res$cfi) & !is.na(res$rmsea) & !is.na(res$srmr)), sum(!is.na(res$p)), sum(!is.na(res$df)))
med_by_review <- tapply(res$n_indices, res$review_id, median)
say("  Median of the per-review medians: %g", median(med_by_review))
write_out(q2, "q2_indices.csv")

# ---- Agreement of printed ratings with the rule as written -------------------------------------
own <- res$rating %in% c("+", "-") & res$rating_scope == "result"
say("")
say("Printed result ratings against the rule as written (results with their own + or - rating, n = %d)",
    sum(own))
agree <- table(printed = res$rating[own], as_written = res$rating_cosmin[own])
for (p in rownames(agree)) say("  printed %s: %s", p, paste(colnames(agree), agree[p, ], sep = " ", collapse = ", "))
plus_no_pass <- own & res$rating == "+" & res$rating_cosmin == "insufficient"
minus_pass <- own & res$rating == "-" & res$rating_cosmin == "sufficient"
say("  Rated + although no reported index meets its threshold: %d results in %d reviews",
    sum(plus_no_pass), length(unique(res$review_id[plus_no_pass])))
say("  Rated - although a reported index meets its threshold: %d results in %d reviews",
    sum(minus_pass), length(unique(res$review_id[minus_pass])))
evaluable <- own & res$rating_cosmin != "indeterminate"
cons <- boot(((res$rating == "+") == (res$rating_cosmin == "sufficient"))[evaluable], res$review_id[evaluable])
say("  Printed rating agrees with the rule as written: %s", cons$text)
dis_rev <- tapply(plus_no_pass | minus_pass, res$review_id, any)[unique(res$review_id[evaluable])]
say("  Reviews with at least one such disagreement: %d of %d reviews with evaluable rated results",
    sum(dis_rev), length(dis_rev))

# ---- Q3 ---------------------------------------------------------------------------------
fails <- function(d) {
  f <- cbind(d$cfi <= 0.95, d$tli <= 0.95, d$rmsea >= 0.06, d$srmr >= 0.08)
  rowSums(f, na.rm = TRUE) > 0
}
res$any_fail <- fails(res)
q3_set <- res$rating == "+" & res$rating_scope == "result" & res$n_indices >= 2
say("")
say("Q3. Results rated sufficient with two or more indices reported: %d (from %d reviews)",
    sum(q3_set), length(unique(res$review_id[q3_set])))
q3 <- boot(res$any_fail[q3_set], res$review_id[q3_set])
say("  At least one reported index fails its own threshold: %s", q3$text)
fail_which <- with(res[q3_set, ], c(CFI = sum(cfi <= .95, na.rm = TRUE), TLI = sum(tli <= .95, na.rm = TRUE),
                                    RMSEA = sum(rmsea >= .06, na.rm = TRUE), SRMR = sum(srmr >= .08, na.rm = TRUE)))
say("  Failing index among these results: %s", paste(names(fail_which), fail_which, sep = " ", collapse = ", "))
say("  Rescued by RMSEA alone (RMSEA < .06 the only passing clause): %d", sum(res$rmsea_rescue[q3_set], na.rm = TRUE))
only_srmr <- with(res, !is.na(srmr) & srmr < .08 & !(cfi > .95 & !is.na(cfi)) & !(tli > .95 & !is.na(tli)) &
                    !(rmsea < .06 & !is.na(rmsea)))
say("  Rescued by SRMR alone: %d", sum(only_srmr[q3_set]))
q3_plus <- res$rating == "+" & res$rating_scope == "result" & res$n_indices >= 1
nopass <- boot((res$rating_cosmin == "insufficient")[q3_plus], res$review_id[q3_plus])
say("  Among all results rated sufficient (one or more indices, n = %d): no reported index passes: %s",
    sum(q3_plus), nopass$text)

# ---- Q4: result level --------------------------------------------------------------------------
say("")
say("Q4a. Re-rating every CFA result under each rule (n = %d results)", nrow(res))
q4 <- do.call(rbind, lapply(RULES, function(r) {
  x <- res[[paste0("rating_", r)]]
  data.frame(rule = r, sufficient = sum(x == "sufficient"), insufficient = sum(x == "insufficient"),
             indeterminate = sum(x == "indeterminate"))
}))
for (i in seq_len(nrow(q4))) say("  %-13s sufficient %3d  insufficient %3d  indeterminate %3d",
                                 q4$rule[i], q4$sufficient[i], q4$insufficient[i], q4$indeterminate[i])
suff <- res$rating_cosmin == "sufficient"
say("  Results sufficient under the rule as written: %d", sum(suff))
q4_change <- do.call(rbind, lapply(REPAIRED, function(r) {
  x <- res[[paste0("rating_", r)]]
  down <- boot((x == "insufficient")[suff], res$review_id[suff])
  lost <- boot((x != "sufficient")[suff], res$review_id[suff])
  up <- sum(res$rating_cosmin == "insufficient" & x == "sufficient")
  say("  %-13s sufficient -> insufficient: %s", r, down$text)
  say("  %-13s sufficient -> not sufficient (insufficient or indeterminate): %s", "", lost$text)
  data.frame(rule = r, n_sufficient_as_written = sum(suff), to_insufficient = sum((x == "insufficient")[suff]),
             to_indeterminate = sum((x == "indeterminate")[suff]), lower = down$lower, upper = down$upper,
             insufficient_to_sufficient = up)
}))
nr <- suff & res$rating_no_rmsea == "insufficient"
say("  no_rmsea downgrades where SRMR was not reported: %d of %d", sum(nr & is.na(res$srmr)), sum(nr))
write_out(q4, "q4_result_ratings.csv")
write_out(q4_change, "q4_result_changes.csv")
by_version <- function(r) {
  x <- res[[paste0("rating_", r)]]
  sapply(1:2, function(v) {
    s <- suff & res$review_id %in% reviews$review_id[reviews$cosmin_version == v]
    sprintf("v%d %d/%d", v, sum((x != "sufficient")[s]), sum(s))
  })
}
for (r in REPAIRED) say("  %-13s not sufficient any more, by COSMIN version: %s", r, paste(by_version(r), collapse = "; "))

# ---- Q4: instrument level ----------------------------------------------------------------------
# An instrument's structural validity under a rule summarises its evaluable results:
# + when at least 75% are sufficient, - when at least 75% are insufficient, +/- otherwise,
# ? when none can be evaluated. Rows are grouped by the instrument the review rated.
pool <- function(x) {
  e <- x[x != "indeterminate"]
  if (!length(e)) return("?")
  s <- mean(e == "sufficient")
  if (s >= .75) "+" else if (s <= .25) "-" else "+/-"
}
res$group <- ifelse(res$instrument_group == "", res$instrument, res$instrument_group)
key <- paste(res$review_id, res$group, sep = "|")
ilev <- data.frame(key = unique(key), stringsAsFactors = FALSE)
ilev$review_id <- sub("\\|.*", "", ilev$key)
ilev$group <- sub("^[^|]*\\|", "", ilev$key)
for (r in RULES) ilev[[paste0("sv_", r)]] <- sapply(ilev$key, function(k) pool(res[[paste0("rating_", r)]][key == k]))
ilev$cosmin_version <- reviews$cosmin_version[match(ilev$review_id, reviews$review_id)]
ilev$scheme <- reviews$recommendation_scheme[match(ilev$review_id, reviews$review_id)]
m <- match(ilev$key, paste(inst$review_id, inst$instrument_group, sep = "|"))
ilev$sv_printed <- inst$sv[m]
ilev$sv_quality <- inst$sv_quality[m]
ilev$ic_printed <- inst$ic[m]
ilev$recommendation <- inst$recommendation[m]
for (v in c("sv_printed", "sv_quality", "ic_printed", "recommendation")) ilev[[v]][is.na(ilev[[v]])] <- ""
say("")
say("Q4b. Instruments (review x instrument): %d, of which %d have a review-level record",
    nrow(ilev), sum(!is.na(m)))
say("  Structural validity + under the rule as written: %d", sum(ilev$sv_cosmin == "+"))
for (r in REPAIRED) {
  x <- ilev[[paste0("sv_", r)]]
  base <- ilev$sv_cosmin == "+"
  say("  %-13s + -> not +: %3d of %d (to -: %d, to +/-: %d, to ?: %d); v1 %d, v2 %d", r,
      sum(base & x != "+"), sum(base), sum(base & x == "-"), sum(base & x == "+/-"), sum(base & x == "?"),
      sum(base & x != "+" & ilev$cosmin_version == 1), sum(base & x != "+" & ilev$cosmin_version == 2))
}

# Internal consistency: COSMIN rates internal consistency only when there is at least low-quality
# evidence of sufficient structural validity (version 1) or of unidimensionality (version 2).
# A printed + or - internal consistency rating becomes ? when structural validity is no longer +.
say("")
say("Q4c. Internal consistency carried forward (instruments with printed IC + or -: %d)",
    sum(ilev$ic_printed %in% c("+", "-")))
ic_tab <- do.call(rbind, lapply(REPAIRED, function(r) {
  x <- ilev[[paste0("sv_", r)]]
  ch <- ilev$ic_printed %in% c("+", "-") & ilev$sv_cosmin == "+" & x != "+"
  say("  %-13s IC -> ?: %d (v1 %d, v2 %d)", r, sum(ch), sum(ch & ilev$cosmin_version == 1),
      sum(ch & ilev$cosmin_version == 2))
  data.frame(rule = r, ic_changed = sum(ch), v1 = sum(ch & ilev$cosmin_version == 1),
             v2 = sum(ch & ilev$cosmin_version == 2))
}))

# Recommendation (version 1, A/B/C): A needs at least low-quality evidence of sufficient internal
# consistency, so A falls to B when internal consistency becomes ?; B becomes C when structural
# validity becomes - with high-quality evidence.
abc <- ilev$cosmin_version == 1 & ilev$scheme == "A/B/C"
say("")
say("Q4d. Recommendations, COSMIN version 1 reviews with A/B/C categories: %d reviews, %d instruments",
    length(unique(ilev$review_id[abc])), sum(abc))
say("  Category printed for %d of these instruments: A %d, B %d, C %d", sum(abc & ilev$recommendation != ""),
    sum(abc & ilev$recommendation == "A"), sum(abc & ilev$recommendation == "B"),
    sum(abc & ilev$recommendation == "C"))
rec_tab <- do.call(rbind, lapply(REPAIRED, function(r) {
  x <- ilev[[paste0("sv_", r)]]
  ic_q <- ilev$ic_printed %in% c("+", "-") & ilev$sv_cosmin == "+" & x != "+"
  a_b <- abc & ilev$recommendation == "A" & ic_q
  b_c <- abc & ilev$recommendation == "B" & ilev$sv_cosmin != "-" & x == "-" & ilev$sv_quality == "high"
  denom <- sum(abc & ilev$recommendation != "")
  w <- wilson_ci(sum(a_b | b_c), denom)
  say("  %-13s A -> B: %d; B -> C: %d; changed %s", r, sum(a_b), sum(b_c),
      sprintf("%d/%d = %s", sum(a_b | b_c), denom, ci_text(w$est, w$lower, w$upper)))
  data.frame(rule = r, a_to_b = sum(a_b), b_to_c = sum(b_c), changed = sum(a_b | b_c), n = denom)
}))
v2 <- ilev$cosmin_version == 2
say("  COSMIN version 2 reviews (%d reviews, %d instruments): structural validity + -> not + under %s",
    length(unique(ilev$review_id[v2])), sum(v2),
    paste(sapply(REPAIRED, function(r) sprintf("%s %d", r, sum(v2 & ilev$sv_cosmin == "+" &
                                                              ilev[[paste0("sv_", r)]] != "+"))), collapse = ", "))
write_out(ilev[, setdiff(names(ilev), "key")], "q4_instruments.csv")
write_out(merge(ic_tab, rec_tab, by = "rule"), "q4_carry_forward.csv")

# ---- Predictions ----------------------------------------------------------------------------
say("")
say("Predictions")
say("  P1 fewer than half of reviews state the criterion unambiguously: %s (%s)",
    if (unamb / nrow(reviews) < .5) "holds" else "does not hold", pct(unamb / nrow(reviews)))
say("  P2 median number of indices per result is 2 or fewer: %s (median %g)",
    if (qs[2] <= 2) "holds" else "does not hold", qs[2])
say("  P3 at least 20%% of sufficient results have a failing index: %s (%s)",
    if (q3$est >= .2) "holds" else "does not hold", pct(q3$est))
max_rec <- max(rec_tab$changed / rec_tab$n)
say("  P4 fewer than 10%% of recommendations change under any repaired rule: %s (largest %s)",
    if (max_rec < .1) "holds" else "does not hold", pct(max_rec))

# ---- Sensitivity analyses ----------------------------------------------------------------------
say("")
say("Sensitivity analyses")
q3_line <- function(d, label) {
  d$any_fail <- fails(d)
  s <- d$rating == "+" & d$rating_scope == "result" & d$n_indices >= 2
  b <- boot(d$any_fail[s], d$review_id[s])
  say("  Q3 %s: %s", label, b$text)
}
change_line <- function(d, label) {
  s <- d$rating_cosmin == "sufficient"
  say("  Q4 %s: sufficient results %d; not sufficient any more: %s", label, sum(s),
      paste(sapply(REPAIRED, function(r) sprintf("%s %d", r, sum(d[[paste0("rating_", r)]][s] != "sufficient"))),
            collapse = ", "))
}
# S1 implausible values set to missing
s1 <- cfa[cfa$level == "study" & !cfa$range, ]
for (k in c("cfi", "tli", "rmsea", "srmr")) s1[[k]][has_flag(s1$value_flags, paste0(k, ":implausible"))] <- NA
s1 <- rate_all(s1)
say("  S1 numerically implausible values treated as not reported (%d values in %d rows)",
    sum(lengths(regmatches(cfa$value_flags, gregexpr("implausible", cfa$value_flags)))), sum(cfa$implausible))
q3_line(s1, "S1"); change_line(s1, "S1")
# S2 rows with printed ranges added: an index passes when both endpoints pass and fails when both fail
endpoints <- function(s) {
  s <- gsub("−", "-", s)
  v <- regmatches(s, gregexpr("[0-9]*\\.?[0-9]+", s))
  t(sapply(v, function(x) if (length(x) == 2) as.numeric(x) else c(NA, NA)))
}
s2 <- cfa[cfa$level == "study", ]
for (k in c("cfi", "tli", "rmsea", "srmr")) {
  rg <- has_flag(s2$value_flags, paste0(k, ":range"))
  if (!any(rg)) next
  e <- endpoints(s2[[paste0(k, "_printed")]][rg])
  lo <- pmin(e[, 1], e[, 2]); hi <- pmax(e[, 1], e[, 2])
  val <- if (k %in% c("cfi", "tli")) ifelse(lo > .95, lo, ifelse(hi <= .95, hi, NA)) else {
    cut <- if (k == "rmsea") .06 else .08
    ifelse(hi < cut, hi, ifelse(lo >= cut, lo, NA))
  }
  s2[[k]][rg] <- val
}
s2 <- rate_all(s2)
say("  S2 study-level rows with printed ranges included (%d rows added)", sum(cfa$level == "study" & cfa$range))
q3_line(s2, "S2"); change_line(s2, "S2")
# S3 comparable incremental indices (IFI, NFI, RNI, RFI) counted with CFI
s3 <- cfa[cfa$level == "study" & !cfa$range, ]
comp <- pmax(s3$ifi, s3$nfi, s3$rni, s3$rfi, na.rm = TRUE)
add <- !is.na(comp) & (is.na(s3$cfi) | comp > s3$cfi)
s3$cfi[add] <- comp[add]
s3 <- rate_all(s3)
say("  S3 IFI, NFI, RNI or RFI > .95 accepted as a comparable measure (%d rows affected)", sum(add))
q3_line(s3, "S3"); change_line(s3, "S3")
# S4 pooled and summary ratings attached to single results included in Q3
s4 <- res
s4$rating_scope[s4$rating_scope == "pooled"] <- "result"
q3_line(s4, "S4 ratings printed at instrument level included")
# S5 an index that was not reported counts as unknown rather than absent
s5 <- rate_all(cfa[cfa$level == "study" & !cfa$range, ], missing = "unknown")
say("  S5 unreported indices treated as unknown: indeterminate results %s",
    paste(sapply(RULES, function(r) sprintf("%s %d", r, sum(s5[[paste0("rating_", r)]] == "indeterminate"))),
          collapse = ", "))
# S6 instrument pooling with a strict consistency rule: + only when every evaluable result is sufficient
strict <- function(x) {
  e <- x[x != "indeterminate"]
  if (!length(e)) return("?")
  if (all(e == "sufficient")) "+" else if (all(e == "insufficient")) "-" else "+/-"
}
sv6 <- sapply(RULES, function(r) sapply(ilev$key, function(k) strict(res[[paste0("rating_", r)]][key == k])))
say("  S6 instrument + only when all evaluable results are sufficient: + as written %d; + -> not +: %s",
    sum(sv6[, "cosmin"] == "+"),
    paste(sapply(REPAIRED, function(r) sprintf("%s %d", r, sum(sv6[, "cosmin"] == "+" & sv6[, r] != "+"))),
          collapse = ", "))

write_out(res[, c("row_id", "review_id", "instrument", "group", "n_indices", "rating", "rating_scope",
                  paste0("rating_", RULES), "rmsea_rescue", "any_fail")], "rerated_results.csv")
writeLines(report, file.path(out_dir, "study3_results.txt"))
cat(report, sep = "\n")
