# Study 1: what the rule allows at population level

Study 1 uses no sampling in part A and one very large sample per condition in part B. It asks how the four fit indices of the COSMIN rule relate to each other, and what "SRMR < 0.08" means in RMSEA units.

## Part A: unweighted least squares (`code/uls_identity.R`)

Fit a CFA by unweighted least squares (ULS) to a population correlation matrix of p items. Write G for the sum of squared residual correlations, B for the sum of squared correlations below the diagonal, D = p(p+1)/2, df for the model degrees of freedom and df_B = p(p-1)/2. Then

- SRMR = sqrt(G / D)
- RMSEA = sqrt(G / df)
- CFI = 1 - G / B
- TLI = 1 - (G / df) / (B / df_B)

The first two expressions are in Shi, Maydeu-Olivares and Rosseel (2020). The four clauses of the rule are therefore four thresholds on the same quantity G, and the rule passes a model when G is below the largest of them:

| Clause | Passes when |
|---|---|
| SRMR < 0.08 | G < 0.0064 D |
| RMSEA < 0.06 | G < 0.0036 df |
| CFI > 0.95 | G < 0.05 B |
| TLI > 0.95 | G < 0.05 B df / df_B |

Because df is smaller than D and than df_B, the RMSEA and TLI clauses are never the most lenient. SRMR = RMSEA sqrt(df / D), so SRMR < 0.08 allows an RMSEA of 0.08 sqrt(D / df).

The script checks these relations in lavaan on 101 misspecified population models (4 to 30 items, one to four factors) and writes Table 1 (`results/table1_rmsea_equivalent.csv`): the RMSEA that SRMR < 0.08 allows for each scale length, and the SRMR cut-off that would hold RMSEA below 0.06.

## Part B: diagonally weighted least squares (`code/dwls_weights.R`)

With ordinal items, DWLS weights each squared residual of a polychoric correlation by w = 1 / (N x its sampling variance). With t = p(p+1)/2 and F the weighted sum of squared residuals,

- chi-square = (N - 1) F and t SRMR^2 = the unweighted sum of squared residuals;
- F = w_eff t SRMR^2, where w_eff is the mean weight, weighted by the squared residuals;
- RMSEA^2 = w_eff (t / df) SRMR^2 - 1 / (N - 1) for the unscaled statistic;
- sqrt(w_min t / df) SRMR <= sqrt(F / df) <= sqrt(w_max t / df) SRMR.

These hold in every sample. With all weights equal to 1 they reduce to the ULS relation of part A. The mean-and-variance scaled statistic divides the chi-square by a scaling factor c, so for the scaled RMSEA that WLSMV reports the weights are in effect w / c.

The script checks the relations and tabulates the weights for the six item formats, two loadings and three scale lengths of Study 2b, with one sample of N = 200,000 per condition. The mean weight w_eff was below 1 in every condition with binary items and smallest with skewed binary items. It rises with the number of categories and with the size of the correlations. The RMSEA that corresponds to a given SRMR therefore depends on the item format as well as on the scale length.

## Files

| File | Contents |
|---|---|
| `code/uls_identity.R` | Part A |
| `code/dwls_weights.R` | Part B |
| `results/uls_identity.txt`, `results/uls_models.csv` | Output of part A: checks, and one row per population model |
| `results/table1_rmsea_equivalent.csv` | Table 1 |
| `results/dwls_weights.txt`, `results/dwls_weights.csv` | Output of part B: checks, and one row per condition |
| `reproduce.R` | Re-runs part A, or both parts with `--dwls`, and compares the values with `results/` |

## Running it

```bash
R CMD INSTALL ../cosminfit
Rscript reproduce.R                 # part A, a few seconds
Rscript reproduce.R --dwls          # parts A and B, about 15 minutes
```

`tests/dwls_identity.R` in the repository root checks the relations of part B on small samples in a few seconds.

## Reference

Shi D, Maydeu-Olivares A, Rosseel Y. Assessing fit in ordinal factor analysis models: SRMR vs. RMSEA. Struct Equ Modeling 2020;27:1-15. https://doi.org/10.1080/10705511.2019.1611434
