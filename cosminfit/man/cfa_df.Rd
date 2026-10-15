\name{cfa_df}
\alias{cfa_df}
\title{Degrees of freedom of a simple-structure CFA}
\usage{
cfa_df(p, m = 1)
}
\arguments{
\item{p}{Number of items (vector allowed).}

\item{m}{Number of factors (vector allowed, recycled with \code{p}).}
}
\value{
Integer-valued numeric vector of degrees of freedom. An error is
raised if any df is negative. A df of 0 gives a warning, because the
SRMR-RMSEA mapping is then undefined.
}
\description{
Degrees of freedom of a CFA model with \code{p} items loading on \code{m}
correlated factors, each item on exactly one factor:
\code{p * (p + 1) / 2 - (2 * p + m * (m - 1) / 2)}.
}
\examples{
cfa_df(4)
cfa_df(10, m = 2)
cfa_df(c(4, 6, 8, 10))
}
