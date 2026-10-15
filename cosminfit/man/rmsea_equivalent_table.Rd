\name{rmsea_equivalent_table}
\alias{rmsea_equivalent_table}
\title{RMSEA equivalent of an SRMR cut-off}
\usage{
rmsea_equivalent_table(
  p = c(4, 5, 6, 8, 10, 15, 20, 30),
  m = 1,
  srmr_cut = 0.08,
  convention = c("lavaan", "offdiag")
)
}
\arguments{
\item{p}{Numbers of items.}

\item{m}{Number of factors (recycled with \code{p}).}

\item{srmr_cut}{SRMR cut-off. Default 0.08 (COSMIN).}

\item{convention}{Passed to \code{\link[=n_moments]{n_moments()}}.}
}
\value{
A data.frame with columns \code{p}, \code{m}, \code{df}, \code{D} and
\code{rmsea_equivalent}.
}
\description{
For each number of items, the RMSEA value that corresponds, at population
level under ULS, to the SRMR cut-off in a one- or multi-factor
simple-structure CFA. See \code{\link[=srmr_to_rmsea]{srmr_to_rmsea()}}.
}
\examples{
rmsea_equivalent_table()
rmsea_equivalent_table(p = c(6, 9, 12), m = 3)
}
