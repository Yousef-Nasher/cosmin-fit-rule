# Study 2: simulation of the rule with ML and WLSMV

How often does the COSMIN rule pass correct and misspecified models, and which clause decides the rating?

## Design

- Scales of 5 items (one factor), 10 items (one or two factors) and 12 items (two factors); loadings of 0.50 or 0.70; two factors correlate 0.30.
- Model conditions: correct; one, two or three residual correlations (0.30, 0.30 and 0.40); two cross-loadings of 0.40; an omitted dimension correlated 0.70. Combinations that do not apply are dropped, which leaves 38 conditions (`code/design.R`).
- N = 100, 200, 500 and 1000, giving 152 cells with 500 replications each.
- Each data set is analysed three ways: ML on the continuous items, ML on the same items cut into five symmetric categories and treated as continuous, and WLSMV on the five-category items. For WLSMV the scaled indices are primary and the robust indices are a sensitivity analysis.
- A misspecification counts as consequential when, at population level, coefficient omega is overstated by 0.05 or more, a factor correlation or a loading is biased by 0.10 or more, or a dimension correlated 0.70 or less is forced into one factor (`code/population.R`). This holds for 21 of the 30 misspecified conditions.

## Files

| File | Contents |
|---|---|
| `code/design.R` | Conditions, population matrices, model syntax, job table and seeds |
| `code/simulate.R`, `code/worker.R` | Data generation and model fitting for one job; the worker loop |
| `code/rerun_job.R` | Re-runs one job from its seed and compares it with the stored file |
| `code/population.R` | Population fit and severity for the 38 conditions |
| `code/rules.R` | The rating rules, shared by the analysis scripts |
| `code/analysis.R`, `code/analysis_extra.R` | Main analysis; effect of the RMSEA clause, results per condition, inadmissible fits kept |
| `code/table_rules.R` | False-fail and false-pass rates of every rule by analysis arm and N |
| `code/hu_bentler_sensitivity.R` | Sensitivity of the Hu-Bentler combination rule to its CFI or TLI cut-off |
| `code/verify_independent.py` | Independent check of the headline numbers in pandas |
| `data/study2_final_chunks.zip` | Stored output: 1,520 job files of 50 replications (228,000 fits) and `population.csv` |
| `results/` | Output of the scripts above |
| `reproduce.R` | Re-runs the analysis from the stored output and compares it with `results/` |

## Running it

```bash
Rscript reproduce.R                 # from this folder
```

To run the simulation itself, start one or more workers from `code/`:

```bash
Rscript worker.R w1                 # writes one CSV per job to code/out/chunks/
```

Each job has the seed `20260927 + 1000 x chunk + cell index`. The column `worker` in the stored files shows where a job ran: `w1` and `w2` used lavaan 0.6-17, and `L1` and `L2` used lavaan 0.6-21.

## Columns of the job files

`cell`, `p`, `m`, `load`, `kind`, `k`, `delta`, `n` identify the cell; `chunk` and `rep` the replication; `arm` the analysis (`ml_cont`, `ml_ord_as_cont`, `wlsmv_ord`). `status` is `ok`, `nonconv`, `degenerate` or `error`, and `admiss` is `TRUE` when no standardised loading or factor correlation reaches 0.999 and no residual variance is negative. The fit indices follow: `chisq`, `df`, `cfi`, `tli`, `rmsea`, `srmr`, then the scaled (`.sc`) and robust (`.rob`) versions for WLSMV. `seed` and `worker` close the row.
