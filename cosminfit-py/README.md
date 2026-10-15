# cosminfit (Python)

Rate confirmatory factor analysis (CFA) results against the COSMIN
structural validity criterion

> CFI or TLI or comparable measure > 0.95 OR RMSEA < 0.06 OR SRMR < 0.08

(Prinsen et al. 2018; Mokkink et al. 2024) and against alternative rules.
The package also gives the population relation between SRMR and RMSEA, so
you can see what RMSEA value the SRMR < 0.08 cut-off corresponds to.

This is a pure-Python port of the R package `cosminfit`. It uses only the
standard library (pandas is optional). Function and argument names are the
same as in R, and the results agree exactly with the R version: the test
suite checks about 150,000 ratings and numbers produced by R.

## Install

```sh
pip install .                            # from this folder
pip install --no-build-isolation .       # offline, uses the installed setuptools
```

## Example

```python
from cosminfit import *
cosmin_rules()                                   # the available rules
rate_cfa(cfi=0.93, rmsea=0.05)                   # 'sufficient': RMSEA alone passes
rate_cfa(cfi=0.93, rmsea=0.05, rule="two_of_three")  # 'insufficient': 1 of 2 reported passes
rate_cfa(cfi=0.93, missing="unknown")            # 'indeterminate'
rmsea_equivalent_table()                         # SRMR .08 as RMSEA, by items
d = [dict(cfi=.97, tli=None, rmsea=.05, srmr=.04, p=10, df=35),
     dict(cfi=.91, tli=.89, rmsea=.05, srmr=.09, p=6, df=9)]
rerate(d)                                        # one rating_<rule> key per rule
wilson_ci(5, 100)                                # Wilson 95% interval
```

`rerate()` also takes a pandas DataFrame and then returns a DataFrame.

## Differences from R

* Scalars in, scalar out: `rate_cfa(cfi=0.97)` returns a string; any list,
  tuple, Series or array returns a list. Tables (`cosmin_rules()`,
  `rmsea_equivalent_table()`, vector `wilson_ci()`) are lists of dicts.
* "Not reported" is `None`, NaN, infinity, an empty string or a string
  that is not a number. A missing numeric result is `float("nan")`.
* `missing`, `rule` and `convention` must be spelled in full (R allows
  partial matching for `missing` and `convention`).
* `cluster_boot_prop()` uses `random.Random(seed)`, so its interval is not
  the same as R's for the same seed (it agrees up to Monte Carlo error).

## Missing indices

Studies often report only some indices. `missing="reported"` (the default)
rates a result using only the indices that were reported.
`missing="unknown"` treats a missing index as possibly passing or failing
and gives a rating only when it does not depend on the missing values. See
`help(rate_cfa)` for the exact logic of each rule.

## Tests

```sh
python -m unittest -v           # or: pytest
```

`tests/test_crosscheck.py` reads the CSV files written by
`../crosscheck/make_vectors.R` and is skipped if they are not there.

## License

MIT
