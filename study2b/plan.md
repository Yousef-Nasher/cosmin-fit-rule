# Study 2b: item format. Analysis plan

Written 29 September 2026, before any Study 2b data were generated. Study 2b is a **secondary study and is not pre-registered**. It was prompted by an exploratory run (about 3,000 WLSMV fits at 5 replications per cell), and the paper will say so. Any change to this plan is recorded with a date and a reason.

## Why

Study 2 uses 5-category symmetric items. In that setting the RMSEA clause of the COSMIN rule almost never changes a rating. The exploratory run suggested this reverses with binary, 3-category or skewed items at small N: SRMR then fails correct models, and RMSEA is what rescues them. Study 2b measures how often this happens, so that the paper can qualify its claim by item format and N rather than state it in general.

## Questions

1. How often does each clause (CFI, TLI, RMSEA, SRMR) fail a **correct** model, by item format, scale length, loading and N?
2. How often is RMSEA the **only** clause that passes a correct model (the RMSEA rescue rate)?
3. How often does the full COSMIN rule pass a model with consequential misspecification?

## Expectations, stated in advance

- E1. With binary or 3-category skewed items at N ≤ 200, SRMR ≥ .08 in more than 20% of correct-model fits.
- E2. In those same cells, the RMSEA rescue rate exceeds 5%.
- E3. With 5-category symmetric items, the RMSEA rescue rate is below 1% in every cell. This would replicate an earlier check of 6,000 fits with 5-category symmetric items and loadings of .70.

We report all three whichever way they come out.

## Design (full factorial, 288 cells)

| Factor | Levels |
|---|---|
| Categories × threshold shape | 2, 3 and 5 categories, each symmetric or skewed (floor effect). Six formats in all. |
| Items (p) | 5, 10, 20 |
| Loading | .45, .70 (all items equal) |
| N | 100, 200, 500 |
| Model | correct; one residual correlation of .30 (items 1 and 2); two factors correlated .70 fitted as one (p ≥ 10 only) |

Marginal category proportions (%):

| Format | Proportions |
|---|---|
| 2 symmetric | 50 / 50 |
| 2 skewed | 88.5 / 11.5 |
| 3 symmetric | 22.7 / 54.7 / 22.7 |
| 3 skewed | 65.5 / 27.8 / 6.7 |
| 5 symmetric | 6.7 / 24.2 / 38.3 / 24.2 / 6.7 |
| 5 skewed | 42.1 / 30.5 / 15.9 / 8.6 / 2.9 |

These are the thresholds from the exploratory run, kept unchanged.

**Replications:** 200 per cell, in 4 jobs of 50. Seed = 20261001 + 1000·chunk + cell index.

## Data and estimation

- Latent responses are drawn from N(0, Σ) and cut at the thresholds. For each data set we fit a one-factor model in two ways:
  - **WLSMV on polychorics.** This is the primary arm. Scaled indices are primary, robust indices are the sensitivity analysis.
  - **ML treating the items as continuous.** This is a secondary arm, because many PROM studies do it.
- Code: `design2b.R`, `simulate2b.R` and `worker2b.R`, with lavaan 0.6-17.
- Unobserved categories are kept as they fall. We do not re-draw data or collapse categories, because a reviewer sees what the data give. The number of affected items is recorded.
- A data set in which any item is constant is marked "degenerate" and not fitted.

## Severity

Severity comes from population values on the latent scale (`population2b.R`), using the same rule as Study 2. A misspecification is consequential if omega is overstated by ≥ .05, any loading is off by ≥ .10, or it is the two-factor condition. This gives:

- the residual correlation is consequential at p = 5 only;
- the two-factor condition is consequential by definition.

Severity never uses fit indices, so the evaluation is not circular.

## Outcomes and analysis

- Primary analysis set: converged and admissible fits (no loading ≥ .999 and no negative residual variance), as in Study 2.
- Sensitivity set: all converged fits.
- We report the rates of non-convergence, inadmissibility and degenerate data per cell.
- For each cell and arm we report:
  - the false-fail rate of each clause and of the COSMIN rule, on correct models;
  - the RMSEA rescue rate: RMSEA < .06 while CFI ≤ .95, TLI ≤ .95 and SRMR ≥ .08;
  - the SRMR-sole rate: SRMR < .08 while the other three fail;
  - the COSMIN pass rate on consequential misspecifications (false pass);
  - how often SRMR > RMSEA.
- All rates come with Wilson 95% intervals. There are no hypothesis tests; this is a description of operating characteristics.
- Main display: one heat map of the RMSEA rescue rate (format × N, faceted by p and loading), plus one table.

## What would count against the paper

If E2 fails, meaning rescue rates stay below 5% even with skewed binary items, then the claim that the RMSEA clause is redundant holds almost everywhere, and the paper says so without qualification. If E3 fails, Study 2's redundancy result is itself fragile, and the discussion must say that.
