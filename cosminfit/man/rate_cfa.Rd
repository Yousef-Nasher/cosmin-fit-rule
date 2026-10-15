\name{rate_cfa}
\alias{rate_cfa}
\title{Rate CFA fit against COSMIN and alternative rules}
\usage{
rate_cfa(
  cfi = NA,
  tli = NA,
  rmsea = NA,
  srmr = NA,
  rule = "cosmin",
  p = NA,
  df = NA,
  missing = c("reported", "unknown")
)
}
\arguments{
\item{cfi, tli, rmsea, srmr}{Fit index values.}

\item{rule}{One rule name from \code{\link[=cosmin_rules]{cosmin_rules()}}. Default \code{"cosmin"}.}

\item{p}{Number of items (only used by \code{srmr_df}).}

\item{df}{Model degrees of freedom (only used by \code{srmr_df}).}

\item{missing}{How to treat indices that were not reported:
\code{"reported"} (default) or \code{"unknown"}.}
}
\value{
Character vector of ratings.
}
\description{
Rates each set of fit indices as \code{"sufficient"}, \code{"insufficient"} or
\code{"indeterminate"} under one rule (see \code{\link[=cosmin_rules]{cosmin_rules()}}). Inputs are
recycled like data.frame columns. Values that are \code{NA}, not finite or
cannot be read as a number count as "not reported".
}
\details{
Each rule is a logic formula over pass indicators (for example
\code{CFI > .95}). An indicator is missing when its index is not reported. For
\code{srmr_df}, the SRMR term is missing also when \code{p} or \code{df} is missing or
\code{df <= 0}.

\strong{\code{missing = "reported"}} (primary analysis: use only the indices the
study reported). The formula is evaluated bottom-up; each part is either
evaluable (TRUE/FALSE) or not:
\itemize{
\item an OR group is evaluable if at least one member is evaluable; it is
TRUE if any evaluable member is TRUE;
\item an AND group is evaluable only if all members are evaluable;
\item "at least k of" is evaluable if at least k members are evaluable; it
is TRUE if at least k evaluable members are TRUE.
}

The rating is sufficient (TRUE), insufficient (FALSE) or indeterminate
(not evaluable). So for \code{cosmin}, \code{no_rmsea} and \code{srmr_df} any reported
passing index gives sufficient; \code{two_of_three} needs two reported indices;
\code{all_three} needs all three; each \code{hu_bentler} conjunct needs both its
indices. \code{prinsen2016_a} needs at least one of CFI/TLI and at least one of
RMSEA/SRMR. \code{prinsen2016_b} is evaluable if SRMR is reported, or if RMSEA
and at least one of CFI/TLI are reported. For example, CFI = .97 and
SRMR = .10 with RMSEA not reported is insufficient: the part
"(CFI or TLI) and RMSEA" cannot be evaluated and is ignored, and the SRMR
part fails. Under \code{missing = "unknown"} the same input is indeterminate.

\strong{\code{missing = "unknown"}} (sensitivity analysis). A missing indicator
may be TRUE or FALSE. All completions are enumerated (at most 2^4); the
rating is sufficient or insufficient only if every completion gives the
same result, otherwise indeterminate. This equals three-valued (Kleene)
logic.

With all indices of a rule reported, both modes give the same rating.
}
\examples{
rate_cfa(cfi = 0.93, rmsea = 0.05)
rate_cfa(cfi = 0.93, rmsea = 0.07, srmr = 0.05, rule = "two_of_three")
rate_cfa(cfi = 0.97, rule = "all_three")
rate_cfa(cfi = 0.93, rule = "cosmin", missing = "unknown")

# Several studies at once
rate_cfa(cfi = c(0.96, 0.90, NA), srmr = c(NA, 0.09, NA))

# From a lavaan fit
if (requireNamespace("lavaan", quietly = TRUE)) {
  fit <- lavaan::cfa("visual =~ x1 + x2 + x3 + x4 + x5 + x6",
                     data = lavaan::HolzingerSwineford1939)
  fm <- lavaan::fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", "df"))
  rate_cfa(cfi = fm["cfi"], tli = fm["tli"], rmsea = fm["rmsea"],
           srmr = fm["srmr"], rule = "srmr_df", p = 6, df = fm["df"])
}
}
