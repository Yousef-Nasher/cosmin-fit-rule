\name{n_moments}
\alias{n_moments}
\title{Number of sample moments used by SRMR}
\usage{
n_moments(p, convention = c("lavaan", "offdiag"))
}
\arguments{
\item{p}{Number of items (vector allowed).}

\item{convention}{\code{"lavaan"} (default) or \code{"offdiag"}.}
}
\value{
Numeric vector \code{D}.
}
\description{
The number of elements \code{D} averaged by SRMR for \code{p} items.
\itemize{
\item \code{"lavaan"}: \code{D = p (p + 1) / 2}, including the diagonal. This is
the convention of lavaan and of Shi, Maydeu-Olivares and Rosseel (2020),
who divide by \eqn{p(p+1)/2} even for correlation structures.
\item \code{"offdiag"}: \code{D = p (p - 1) / 2}, correlations below the diagonal only
(an alternative, for sensitivity checks).
}
}
\examples{
n_moments(10)
n_moments(10, convention = "offdiag")
}
