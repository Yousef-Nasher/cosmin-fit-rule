\name{cosmin_rules}
\alias{cosmin_rules}
\title{Available rating rules}
\usage{
cosmin_rules()
}
\value{
A data.frame with columns \code{rule}, \code{definition}, \code{indices} (the
indices the rule uses) and \code{needs_p_df} (whether \code{p} and \code{df} are
needed).
}
\description{
Lists the rules that \code{\link[=rate_cfa]{rate_cfa()}} and \code{\link[=rerate]{rerate()}} can apply. All
inequalities are as written (strict unless shown otherwise).
\itemize{
\item \code{cosmin}: CFI > .95 OR TLI > .95 OR RMSEA < .06 OR SRMR < .08
(Prinsen et al. 2018; Mokkink et al. 2024).
\item \code{no_rmsea}: CFI > .95 OR TLI > .95 OR SRMR < .08.
\item \code{srmr_df}: CFI > .95 OR TLI > .95 OR RMSEA < .06 OR
SRMR < .06 * sqrt(df / D), with D = p(p+1)/2. Needs \code{p} and \code{df}.
\item \code{two_of_three}: at least two of CFI > .95, RMSEA < .06, SRMR < .08
(TLI not used).
\item \code{all_three}: CFI > .95 AND RMSEA < .06 AND SRMR < .08.
\item \code{hu_bentler}: (CFI >= .95 AND SRMR <= .09) OR
(RMSEA <= .06 AND SRMR <= .09).
\item \code{prinsen2016_a}, \code{prinsen2016_b}: two readings of the 2016
COSMIN-COMET guideline (Prinsen et al. 2016, Trials 17:449, Table 4),
which printed "CFI or TLI or comparable measure > 0.95 AND RMSEA < 0.06
OR SRMR < 0.08". Reading a is (CFI > .95 OR TLI > .95) AND
(RMSEA < .06 OR SRMR < .08); reading b is
((CFI > .95 OR TLI > .95) AND RMSEA < .06) OR SRMR < .08.
}
}
\examples{
cosmin_rules()
}
