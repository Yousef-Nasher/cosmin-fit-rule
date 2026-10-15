# cosmin-fit-rule

Code, data and software for a study of the COSMIN criterion for structural validity.

COSMIN rates the structural validity of an outcome measurement instrument as sufficient when a confirmatory factor analysis (CFA) gives "CFI or TLI or comparable measure > 0.95 OR RMSEA < 0.06 OR SRMR < 0.08" (Prinsen et al. 2018; Mokkink et al. 2024). This repository holds everything needed to check what that rule allows, how it behaves in simulated data, and how published reviews apply it.

## What is here

| Folder | Contents |
|---|---|
| `study1/` | Population algebra. Part A: under unweighted least squares, SRMR, RMSEA, CFI and TLI are four functions of the same sum of squared residuals, so the rule is one threshold on that sum. Part B: the same relation under diagonally weighted least squares, where the weights depend on the item format. |
| `study2/` | Main simulation: 152 cells with 500 replications each, analysed with ML and WLSMV in lavaan. Design, simulation code, stored output and analysis. |
| `study2b/` | Second simulation on item format: 288 cells with 200 replications each (binary, three-category and five-category items, symmetric or skewed). Not registered; the plan was written before the data were generated. |
| `study3/` | Meta-research on a stratified random sample of published COSMIN reviews: sampling frame, random draw, full-text status of every record, extracted data, analysis and independent check. Registered on OSF: <https://osf.io/46zdg/overview>. |
| `cosminfit/` | R package that rates CFA results under the COSMIN rule and seven other rules, gives the SRMR and RMSEA mapping, and computes Wilson and cluster-bootstrap intervals. |
| `cosminfit-py/` | The same package in Python, using the standard library only. |
| `crosscheck/` | Test vectors written by the R package and read by the Python tests, so that the two packages are held to the same answers. |
| `figures/` | Script and output for the figures of Studies 2, 2b and 3. |
| `tests/` | Checks on the repository as a whole: checksums of the stored data, consistency between files, and agreement between the R results and the independent Python results. |

Each study folder has its own `README.md`.

## Authors

Yousef M. Nasher<sup>1</sup>, Marwan Alqumbaey<sup>1</sup>, Fareed AL-Mekhlafi<sup>1</sup>, Suhaib Nasr<sup>1</sup>, Hael Saeed Hassan Ali<sup>1</sup>, Abdulaleem Badi<sup>2</sup>, Zaid Al-hayas<sup>1</sup>, Ammar AL-Hadidah<sup>1</sup>

<sup>1</sup> Faculty of Medicine and Health Sciences, Sana'a University, Sana'a, Yemen

<sup>2</sup> Department of Diagnostic Radiography, University of Science and Technology, Sana'a, Yemen

Contact: Yousef M. Nasher, yousefnasher76@gmail.com

## Reproducing the results

Install the R package first, then run the `reproduce.R` script of each study from inside its folder. Each script re-runs the analysis and compares every result with the stored copy. It prints one line per file and stops with an error status if anything differs.

```bash
R CMD INSTALL cosminfit

cd study1  && Rscript reproduce.R && cd ..    # needs lavaan; add --dwls for part B (about 15 minutes)
cd study2  && Rscript reproduce.R && cd ..    # base R; about one minute
cd study2b && Rscript reproduce.R && cd ..    # base R; about one minute
cd study3  && Rscript reproduce.R && cd ..    # base R and cosminfit
```

Studies 2, 2b and 3 also have checks written in Python that share no code with the R analysis. `reproduce.R` runs them when Python with pandas and numpy is available and says so when it is not.

The simulations themselves take many CPU hours, so the repository stores their output (`study2/data/study2_final_chunks.zip` and `study2b/data/study2b_chunks.zip`). Every simulation job has its own seed, and `rerun_job.R` (Study 2) or `rerun_job2b.R` (Study 2b) re-runs a single job and compares it with the stored file.

To redraw the figures, run `python3 figures/make_figures.py` from the repository root.

## Software

The stored results were produced with R 4.3.3 and lavaan 0.6-17. In Study 2, 312 of the 1,520 simulation jobs (workers `L1` and `L2` in the data) were run with lavaan 0.6-21. The Python checks and figures use Python 3 with the package versions in `requirements.txt`.

## Tests

The workflow in `.github/workflows/checks.yml` runs on every push. It builds and checks the R package, runs the Python tests against the cross-check vectors, runs the `reproduce.R` script of each study, and runs the repository tests in `tests/`. `SHA256SUMS` lists the checksums of the stored data files.

## Using the packages

```r
# R
remotes::install_github("Yousef-Nasher/cosmin-fit-rule", subdir = "cosminfit")

library(cosminfit)
rate_cfa(cfi = 0.93, rmsea = 0.05)                      # "sufficient": RMSEA alone passes
rate_cfa(cfi = 0.93, rmsea = 0.05, rule = "no_rmsea")   # "insufficient"
rmsea_equivalent_table()                                # SRMR 0.08 expressed as RMSEA
```

```bash
# Python
pip install "git+https://github.com/Yousef-Nasher/cosmin-fit-rule.git#subdirectory=cosminfit-py"
```

## Licence

Code is released under the MIT licence (see `LICENSE`). Data and text are released under CC BY 4.0.

## Citation

See `CITATION.cff`.

## References

Mokkink LB, Elsman EBM, Terwee CB. COSMIN guideline for systematic reviews of patient-reported outcome measures version 2.0. Qual Life Res 2024;33:2929-39. https://doi.org/10.1007/s11136-024-03761-6

Prinsen CAC, Mokkink LB, Bouter LM, Alonso J, Patrick DL, de Vet HCW, et al. COSMIN guideline for systematic reviews of patient-reported outcome measures. Qual Life Res 2018;27:1147-57. https://doi.org/10.1007/s11136-018-1798-3
