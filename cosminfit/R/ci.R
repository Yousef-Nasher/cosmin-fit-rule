#' Wilson score confidence interval for a proportion
#'
#' @param x Number of successes (vector allowed).
#' @param n Number of trials (vector allowed).
#' @param conf Confidence level. Default 0.95.
#' @return A data.frame with columns `est` (`x / n`), `lower` and `upper`.
#' @examples
#' wilson_ci(5, 100)
#' wilson_ci(c(0, 10, 50), 50)
#' @export
wilson_ci <- function(x, n, conf = 0.95) {
  if (any(n <= 0) || any(x < 0) || any(x > n)) {
    stop("Need n > 0 and 0 <= x <= n.", call. = FALSE)
  }
  check_conf(conf)
  z <- qnorm(1 - (1 - conf) / 2)
  centre <- (x + z^2 / 2) / (n + z^2)
  half <- z * sqrt(x * (n - x) / n + z^2 / 4) / (n + z^2)
  data.frame(est = x / n, lower = centre - half, upper = centre + half)
}

#' Cluster bootstrap confidence interval for a proportion
#'
#' Percentile bootstrap interval for the proportion of `TRUE` values in `x`,
#' resampling whole clusters (for example, studies that report several CFA
#' models) with replacement.
#'
#' @param x Logical (or 0/1) vector. `NA` values are dropped.
#' @param cluster Cluster identifier, same length as `x`.
#' @param B Number of bootstrap samples.
#' @param seed Optional random seed. The global random number state is
#'   restored on exit.
#' @param conf Confidence level. Default 0.95.
#' @return A list with `est` (`mean(x)`), `lower`, `upper`, `B`, `n`
#'   (number of non-missing values) and `n_clusters`.
#' @examples
#' x <- c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)
#' study <- c(1, 1, 2, 3, 3, 3, 4, 5)
#' cluster_boot_prop(x, study, B = 500, seed = 1)
#' @export
cluster_boot_prop <- function(x, cluster, B = 2000, seed = NULL, conf = 0.95) {
  if (length(x) != length(cluster)) {
    stop("`x` and `cluster` must have the same length.", call. = FALSE)
  }
  check_conf(conf)
  keep <- !is.na(x) & !is.na(cluster)
  x <- as.numeric(x[keep])
  cluster <- cluster[keep]
  if (length(x) == 0) stop("No non-missing values in `x`.", call. = FALSE)

  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    }
    on.exit({
      if (is.null(old)) rm(".Random.seed", envir = globalenv())
      else assign(".Random.seed", old, envir = globalenv())
    })
    set.seed(seed)
  }

  successes <- tapply(x, cluster, sum)
  sizes <- tapply(x, cluster, length)
  k <- length(sizes)
  boot <- replicate(B, {
    pick <- sample.int(k, k, replace = TRUE)
    sum(successes[pick]) / sum(sizes[pick])
  })
  alpha <- (1 - conf) / 2
  ci <- quantile(boot, c(alpha, 1 - alpha), names = FALSE)
  list(est = mean(x), lower = ci[1], upper = ci[2], B = B,
       n = length(x), n_clusters = k)
}

check_conf <- function(conf) {
  if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) {
    stop("`conf` must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}
