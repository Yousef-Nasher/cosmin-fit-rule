# cosminfit

Rate confirmatory factor analysis (CFA) results against the COSMIN
structural validity criterion

> CFI or TLI or comparable measure > 0.95 OR RMSEA < 0.06 OR SRMR < 0.08

(Prinsen et al. 2018; Mokkink et al. 2024) and against alternative rules.
The package also gives the population relation between SRMR and RMSEA, so
you can see what RMSEA value the SRMR < 0.08 cut-off corresponds to.

## Install

```r
# install.packages("remotes")
remotes::install_github("Yousef-Nasher/cosmin-fit-rule", subdir = "cosminfit")
```

## Example

```r
library(cosminfit)
cosmin_rules()                                  # the available rules
rate_cfa(cfi = 0.93, rmsea = 0.05)              # "sufficient": RMSEA alone passes
rate_cfa(cfi = 0.93, rmsea = 0.05, rule = "two_of_three")   # "insufficient"
rate_cfa(cfi = 0.93, rmsea = 0.05, rule = "two_of_three",
         missing = "unknown")                                # "indeterminate"
rate_cfa(cfi = 0.93, missing = "unknown")       # "indeterminate"
rmsea_equivalent_table()                        # SRMR .08 as RMSEA, by items
d <- data.frame(cfi = c(.97, .91), tli = c(NA, .89), rmsea = c(.05, .05),
                srmr = c(.04, .09), p = c(10, 6), df = c(35, 9))
rerate(d)                                       # one column per rule
wilson_ci(5, 100)                               # Wilson 95% interval
```

## Missing indices

Studies often report only some indices. `missing = "reported"` (the
default) rates a result using only the indices that were reported.
`missing = "unknown"` treats a missing index as possibly passing or failing
and gives a rating only when it does not depend on the missing values. See
`?rate_cfa` for the exact logic of each rule.

## Citation

Nasher YM, Alqumbaey M, AL-Mekhlafi F, Nasr S, Ali HSH, Badi A, Al-hayas Z, AL-Hadidah A. One passing index is enough: the COSMIN criterion for structural validity in theory, simulation and practice. 2026 (manuscript).

Software archive: https://doi.org/10.5281/zenodo.23405152. See
`citation("cosminfit")`.

## License

MIT
