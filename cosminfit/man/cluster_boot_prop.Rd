\name{cluster_boot_prop}
\alias{cluster_boot_prop}
\title{Cluster bootstrap confidence interval for a proportion}
\usage{
cluster_boot_prop(x, cluster, B = 2000, seed = NULL, conf = 0.95)
}
\arguments{
\item{x}{Logical (or 0/1) vector. \code{NA} values are dropped.}

\item{cluster}{Cluster identifier, same length as \code{x}.}

\item{B}{Number of bootstrap samples.}

\item{seed}{Optional random seed. The global random number state is
restored on exit.}

\item{conf}{Confidence level. Default 0.95.}
}
\value{
A list with \code{est} (\code{mean(x)}), \code{lower}, \code{upper}, \code{B}, \code{n}
(number of non-missing values) and \code{n_clusters}.
}
\description{
Percentile bootstrap interval for the proportion of \code{TRUE} values in \code{x},
resampling whole clusters (for example, studies that report several CFA
models) with replacement.
}
\examples{
x <- c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
study <- c(1, 1, 2, 3, 3, 3, 4, 5)
cluster_boot_prop(x, study, B = 500, seed = 1)
}
