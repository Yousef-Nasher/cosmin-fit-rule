# Study 3 codebook

An empty cell means the item was not reported. We imputed nothing and did not carry values over from another subscale, model or instrument. Review identifiers (`F0001` to `F1177`) are the record numbers of the sampling frame.

## Coding rules

### Criterion wording

We copied the sentence in which a review states the structural validity criterion into `criterion_verbatim` and coded it in `wording_code` as follows.

| Code | Label | Meaning |
|---|---|---|
| 1 | disjunctive | Any one index passing is sufficient ("or" between all indices). |
| 2 | conjunctive | All reported indices must pass ("and", or a stated requirement that indices agree). |
| 3 | mixed | Both connectives in a way that cannot be resolved. A list of thresholds with no connective at all ("CFI > 0.95, RMSEA < 0.06") is coded here. |
| 4 | availability | Rated on whatever the primary study reported, with no stated combination rule. |
| 5 | not stated | No criterion given, or only a citation to COSMIN. |

When a review states the criterion in one place and applies it differently elsewhere, we coded the stated wording and described the discrepancy in `wording_note` or `extraction_note`.

### Fit index values

Values are recorded as printed (`*_printed`). A range is kept as printed and flagged. A value given only as "acceptable" or "good", or as an inequality at the threshold such as "RMSEA < 0.06", counts as not reported. A bound such as "RMSEA < 0.001" is kept and flagged.

### Ratings

Ratings are recorded as printed (`rating_printed`, `sv_printed`, `ic_printed`, `recommendation_printed`) and in a normalised form: `+` sufficient, `-` insufficient, `+/-` inconsistent and `?` indeterminate. We read "supported" and "not supported" as `+` and `-`, and `0` as `?`.

### COSMIN version

The version is 2 when the review follows the COSMIN guideline version 2.0 (Mokkink et al. 2024), and 1 when it follows the 2018 guideline or earlier COSMIN tools.

## data/study3_frame.csv

| Variable | Meaning |
|---|---|
| `frame_id` | Record number |
| `rayyan_key` | Key of the record in the screening export |
| `title`, `authors`, `journal`, `doi`, `url`, `pubmed_id` | Bibliographic data as exported |
| `year_original` | Publication year in the export |
| `year_eff` | Year used for stratification: `year_original`, or the publication year looked up from the DOI or URL when the export had none |
| `year_source` | `rayyan_year` or `recovered_doi_url` |

`data/study3_frame_outside_years.csv` lists the records removed from the frame because `year_eff` was before 2018.

## data/study3_search_counts.csv

`source`, `stage` and `n` for every step from the search to the sampling frame.

## data/study3_fulltext_screening.csv

| Variable | Meaning |
|---|---|
| `frame_id` | Record number |
| `year_stratum` | Publication year |
| `draw_list` | `drawn` (one of the 100) or `replacement queue` |
| `rank` | Draw rank (1 to 100) for drawn records; rank within the stratum for the queue |
| `fulltext_status` | `included`, `excluded`, `full text not obtained` or `not reached` (the quota of the stratum was filled before this rank) |
| `reason` | Reason for exclusion |

## data/study3_reviews.csv

| Variable | Meaning |
|---|---|
| `review_id` | Record number |
| `first_author`, `year`, `title`, `journal`, `doi` | Bibliographic data |
| `year_stratum`, `entered_as` | Stratum, and whether the review was drawn or entered as a replacement |
| `cosmin_version` | 1 or 2 (see above) |
| `cosmin_guideline` | COSMIN publications the review cites for its methods |
| `recommendation_scheme` | `A/B/C` (version 1 categories), `v2`, `narrative` or `none` |
| `criterion_verbatim`, `criterion_location` | The criterion sentence and where it was found |
| `wording_code`, `wording`, `wording_note` | Coding of the criterion (see above) |
| `n_cfa_rows`, `n_instruments` | Rows in the two files below |
| `extraction_note` | Notes on how the review applies the rule |

## data/study3_cfa_results.csv

| Variable | Meaning |
|---|---|
| `review_id`, `row_id` | Review and row |
| `instrument` | Instrument, version or subscale as named in the review |
| `instrument_group` | The instrument as the review rates it at instrument level (empty when the review gives no instrument-level rating) |
| `model`, `study` | Model and primary study, as described |
| `n`, `items` | Sample size and number of items, if reported |
| `estimator`, `ordinal` | Estimator, and whether items were treated as ordinal, if stated (`NS` not stated) |
| `level` | `study` (one CFA result), `pooled` (a result the review pooled over studies) or `summary` (a review-level summary) |
| `cfi_printed` ... `srmr_printed` | Index as printed |
| `cfi`, `tli`, `rmsea`, `srmr` | Index as a number; empty when not reported or printed as a range |
| `ifi`, `nfi`, `rni`, `rfi` | Other incremental indices, if reported |
| `df` | Model degrees of freedom, if reported |
| `other_indices` | Other indices and remarks as printed |
| `rating_printed`, `rating` | Structural validity rating the review gave this result |
| `rating_scope` | `result` when the rating belongs to this result; `pooled` when the review printed an instrument-level rating next to it |
| `value_flags` | `range`, `bound` or `implausible` (CFI or TLI outside 0.5 to 1.2; RMSEA or SRMR above 0.3), per index |
| `source` | Where in the review or its supplement the values were found |

## data/study3_instruments.csv

| Variable | Meaning |
|---|---|
| `review_id`, `instrument_group` | Review and instrument |
| `sv_printed`, `sv`, `sv_quality` | Overall structural validity rating and its quality of evidence |
| `ic_printed`, `ic`, `ic_quality` | Overall internal consistency rating and its quality of evidence. Category A requires sufficient internal consistency, so `ic` is `+` for instruments in category A when the review printed no separate rating |
| `recommendation_printed`, `recommendation` | Recommendation as printed, and the version 1 category (A, B or C) |
