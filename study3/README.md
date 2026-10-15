# Study 3: how published COSMIN reviews apply the rule

This study asks how systematic reviews that use COSMIN state and apply the structural validity criterion, and what would change if the criterion were repaired. It was registered on OSF before any review was read: <https://osf.io/46zdg/overview>.

## Questions

1. Q1. How do reviews state the criterion? We coded it as disjunctive, conjunctive, mixed, an availability rule, or not stated.
2. Q2. How many of the four indices (CFI, TLI, RMSEA, SRMR) does a review report for each CFA result?
3. Q3. Among CFA results rated sufficient, how often does at least one reported index fail its own threshold?
4. Q4. If the rule is repaired, how many structural validity ratings, internal consistency ratings and recommendations change?

## Sample

We searched PubMed with `(COSMIN[tiab]) AND (systematic review[pt] OR systematic review[tiab]) AND ("2018"[dp] : "2026"[dp])` and exported the COSMIN database of systematic reviews. Both were run on 7 October 2026, and the counts for each step are in `data/study3_search_counts.csv`.

The sampling frame (`data/study3_frame.csv`) holds the reviews published from 2018 to 2026 that passed title and abstract screening. Records removed because their publication year was before 2018 are listed in `data/study3_frame_outside_years.csv`.

We drew 100 reviews, stratified by publication year with proportional allocation (largest remainder) and seed 545749, using `code/draw_sample.py`. The records that were not drawn form a ranked replacement queue within each stratum.

A review was included if it applied COSMIN criteria and printed at least one numeric CFA fit index in the main text or in a supplementary file. A review that was ineligible, or whose full text could not be obtained, was replaced by the next record in the queue of the same stratum. The draw order, full-text status and reason for exclusion of every record are in `data/study3_fulltext_screening.csv`.

## Rules compared

| Rule | Sufficient when |
|---|---|
| As written (`cosmin`) | CFI > .95 or TLI > .95 or RMSEA < .06 or SRMR < .08 |
| RMSEA clause dropped (`no_rmsea`) | CFI > .95 or TLI > .95 or SRMR < .08 |
| SRMR scaled by df (`srmr_df`) | CFI > .95 or TLI > .95 or RMSEA < .06 or SRMR < .06 sqrt(df / D), with D = p(p + 1)/2 |
| Two of three (`two_of_three`) | at least two of CFI > .95, RMSEA < .06, SRMR < .08 |

Each result is re-rated with the indices the review reported. If a result has too few indices for a rule, it is indeterminate under that rule, not insufficient. For `srmr_df` we use the reported df. When only the number of items p is known, we use the one-factor value p(p - 3)/2, which is the largest df any model on those items can have and so gives the most lenient cut-off. When SRMR is below .08 but neither p nor df is known, the scaled clause cannot be evaluated, and a result that would otherwise fail is indeterminate.

An instrument's structural validity under a rule is + when at least 75% of its evaluable results are sufficient, - when at most 25% are, +/- otherwise, and ? when no result can be evaluated. A printed internal consistency rating of + or - becomes ? when structural validity is no longer +. In reviews that use the COSMIN version 1 recommendation categories, A becomes B when internal consistency becomes ?, and B becomes C when structural validity becomes - with high-quality evidence. Reviews that follow COSMIN version 2.0 are reported separately.

## Files

| File | Contents |
|---|---|
| `codebook.md` | Variables of every data file and the coding rules |
| `data/study3_search_counts.csv` | Records found, removed and screened before the frame |
| `data/study3_frame.csv` | The sampling frame (1,120 reviews) |
| `data/study3_frame_outside_years.csv` | Records removed from the frame because of a pre-2018 year |
| `data/study3_fulltext_screening.csv` | Draw order and full-text status of every record in the frame |
| `data/study3_reviews.csv` | One row per included review: COSMIN version, verbatim criterion and its coding |
| `data/study3_cfa_results.csv` | One row per CFA result: indices as printed and as numbers, and the rating printed by the review |
| `data/study3_instruments.csv` | One row per instrument rated by a review: structural validity, internal consistency, recommendation |
| `code/draw_sample.py` | The stratified draw |
| `code/analysis.R` | Q1 to Q4, predictions and sensitivity analyses |
| `code/exploratory_q1_period.R` | Exploratory analysis that was not registered: Q1 by publication period (2018 to 2024 and 2025 to 2026) |
| `code/verify_independent.py` | Independent check of the main numbers, written in pandas |
| `results/` | Output of the scripts above |
| `reproduce.R` | Re-runs the draw and the analysis and compares the output with `results/` |

## Running it

From this folder, with the cosminfit package installed:

```bash
Rscript reproduce.R
```

Proportions of reviews are given with Wilson intervals. Proportions of CFA results are given with percentile intervals from a bootstrap that resamples whole reviews (2,000 samples, seed 545749).
