\name{wilson_ci}
\alias{wilson_ci}
\title{Wilson score confidence interval for a proportion}
\usage{
wilson_ci(x, n, conf = 0.95)
}
\arguments{
\item{x}{Number of successes (vector allowed).}

\item{n}{Number of trials (vector allowed).}

\item{conf}{Confidence level. Default 0.95.}
}
\value{
A data.frame with columns \code{est} (\code{x / n}), \code{lower} and \code{upper}.
}
\description{
Wilson score confidence interval for a proportion
}
\examples{
wilson_ci(5, 100)
wilson_ci(c(0, 10, 50), 50)
}
