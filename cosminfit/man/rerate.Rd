\name{rerate}
\alias{rerate}
\title{Re-rate a data set of CFA results under several rules}
\usage{
rerate(
  data,
  rules = c("cosmin", "no_rmsea", "srmr_df", "two_of_three", "all_three"),
  missing = "reported",
  cols = list(cfi = "cfi", tli = "tli", rmsea = "rmsea", srmr = "srmr", p =
    "p", df = "df")
)
}
\arguments{
\item{data}{A data.frame with one row per CFA result.}

\item{rules}{Rule names from \code{\link[=cosmin_rules]{cosmin_rules()}}.}

\item{missing}{\code{"reported"} (default) or \code{"unknown"}; see \code{\link[=rate_cfa]{rate_cfa()}}.}

\item{cols}{Named list giving the column names for \code{cfi}, \code{tli},
\code{rmsea}, \code{srmr}, \code{p} and \code{df}. Partial lists are filled in with the
defaults. A column that is absent from \code{data} counts as not reported.}
}
\value{
\code{data} with one column \verb{rating_<rule>} per rule, plus
\code{n_indices} (how many of CFI, TLI, RMSEA and SRMR were reported) and
\code{rmsea_rescue} (TRUE when RMSEA < .06 is the only one of the four
COSMIN clauses that passes among reported indices; NA when RMSEA was
not reported).
}
\description{
Applies \code{\link[=rate_cfa]{rate_cfa()}} with several rules to each row of a data.frame.
}
\examples{
d <- data.frame(cfi = c(0.97, 0.91, NA), tli = c(NA, 0.89, NA),
                rmsea = c(0.05, 0.05, NA), srmr = c(0.04, 0.09, 0.07),
                p = c(10, 6, 8), df = c(35, 9, 20))
rerate(d)
rerate(d, rules = c("cosmin", "hu_bentler"), missing = "unknown")
}
