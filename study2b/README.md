# Study 2b: item format

Study 2 used five-category symmetric items. Study 2b asks how the clauses of the COSMIN rule behave with binary, three-category and five-category items, with symmetric or skewed thresholds.

This study was not registered. Its plan (`plan.md`) was written on 29 September 2026, before any Study 2b data were generated, and states three expectations. The analysis follows the plan without change.

## Design

- Six item formats: 2, 3 or 5 categories, each symmetric or skewed towards the lowest category (thresholds in `code/design2b.R`).
- 5, 10 or 20 items; loadings of 0.45 or 0.70; N = 100, 200 or 500.
- Three model conditions, all fitted as one factor: correct; one residual correlation of 0.30; two factors correlated 0.70 (10 and 20 items only).
- 288 cells with 200 replications each. Every data set is analysed with WLSMV and with ML on the items treated as continuous.

## Files

| File | Contents |
|---|---|
| `plan.md` | Analysis plan and expectations |
| `code/design2b.R` | Thresholds, conditions, population matrices, job table and seeds |
| `code/simulate2b.R`, `code/worker2b.R` | Data generation and model fitting for one job; the worker loop |
| `code/rerun_job2b.R` | Re-runs one job from its seed and compares it with the stored file |
| `code/population2b.R` | Population fit and severity on the latent-response scale |
| `code/analysis2b.R` | Analysis as fixed in the plan |
| `code/verify_independent_2b.py` | Independent check of the expectation counts in pandas |
| `data/study2b_chunks.zip` | Stored output: 1,152 job files of 50 replications (115,200 fits) and `population2b.csv` |
| `results/` | Output of the analysis and of the independent check |
| `reproduce.R` | Re-runs the analysis from the stored output and compares it with `results/` |

## Running it

```bash
Rscript reproduce.R                 # from this folder
```

Each job has the seed `20261001 + 1000 x chunk + cell index`. All jobs ran with lavaan 0.6-17; the column `lavaan` in the stored files records it. The job files have the columns of Study 2, with `ct` for the item format and `n_empty` for the number of items with an unobserved category.
