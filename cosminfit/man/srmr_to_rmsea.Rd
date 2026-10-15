\name{srmr_to_rmsea}
\alias{srmr_to_rmsea}
\alias{rmsea_to_srmr}
\title{Population relation between SRMR and RMSEA}
\usage{
srmr_to_rmsea(srmr, p, df, convention = c("lavaan", "offdiag"))

rmsea_to_srmr(rmsea, p, df, convention = c("lavaan", "offdiag"))
}
\arguments{
\item{srmr, rmsea}{Index values (vectors allowed).}

\item{p}{Number of items.}

\item{df}{Model degrees of freedom (must be positive), e.g. from
\code{\link[=cfa_df]{cfa_df()}}.}

\item{convention}{Passed to \code{\link[=n_moments]{n_moments()}}: \code{"lavaan"} (default) or
\code{"offdiag"}.}
}
\value{
Numeric vector. \code{NA} (with a warning) where \code{df <= 0}.
}
\description{
Let \code{G} be the sum of squared residual correlations. Then
\code{SRMR = sqrt(G / D)} and, at population level under unweighted least
squares (ULS), \code{RMSEA = sqrt(G / df)}. Hence
\code{RMSEA = SRMR * sqrt(D / df)} and \code{SRMR = RMSEA * sqrt(df / D)}, where
\code{D} is given by \code{\link[=n_moments]{n_moments()}}.
}
\details{
This holds for population values under ULS. It does not hold exactly for
sample values or for other estimators (such as ML), so do not use it to
convert one reported index into another.
}
\examples{
srmr_to_rmsea(0.08, p = 4, df = 2)
rmsea_to_srmr(0.06, p = 10, df = cfa_df(10))
}
