"""Confidence intervals for proportions (port of ``R/ci.R``)."""

from __future__ import annotations

import math
import numbers
import random
from statistics import NormalDist
from typing import Any, Dict, Hashable, List, Optional, Sequence, Union

from ._rules import _as_list, _is_vector

__all__ = ["wilson_ci", "cluster_boot_prop"]


def _check_conf(conf: Any) -> float:
    ok = (not isinstance(conf, bool) and isinstance(conf, numbers.Real)
          and 0 < conf < 1)
    if not ok:
        raise ValueError("`conf` must be a single number between 0 and 1.")
    return float(conf)


def _qnorm_upper(conf: float) -> float:
    """``qnorm(1 - (1 - conf) / 2)`` (Wichura's AS241, as in R)."""
    return NormalDist().inv_cdf(1 - (1 - conf) / 2)


def _num(v: Any) -> float:
    """A count for :func:`wilson_ci`; anything non-numeric or missing errors."""
    if v is None or isinstance(v, (str, bytes, bytearray)):
        raise ValueError("Need n > 0 and 0 <= x <= n.")
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise ValueError("Need n > 0 and 0 <= x <= n.") from None
    if math.isnan(f):
        raise ValueError("Need n > 0 and 0 <= x <= n.")
    return f


def wilson_ci(x: Any, n: Any, conf: float = 0.95
              ) -> Union[Dict[str, float], List[Dict[str, float]]]:
    """Wilson score confidence interval for a proportion.

    Parameters
    ----------
    x
        Number of successes (scalar or sequence).
    n
        Number of trials (scalar or sequence, recycled with ``x``).
    conf
        Confidence level. Default 0.95.

    Returns
    -------
    dict or list of dict
        Keys ``est`` (``x / n``), ``lower`` and ``upper``. A single dict
        when ``x`` and ``n`` are scalars, else one dict per element (the R
        data.frame, row by row).

    Raises
    ------
    ValueError
        Unless ``n > 0`` and ``0 <= x <= n`` everywhere, or if ``conf`` is
        not a single number strictly between 0 and 1.

    Examples
    --------
    >>> w = wilson_ci(5, 100)
    >>> round(w["lower"], 4), round(w["upper"], 4)
    (0.0215, 0.1118)
    """
    scalar = not (_is_vector(x) or _is_vector(n))
    xs = [_num(v) for v in _as_list(x)]
    ns = [_num(v) for v in _as_list(n)]
    if not xs or not ns:
        return []
    k = max(len(xs), len(ns))
    xs = [xs[i % len(xs)] for i in range(k)]
    ns = [ns[i % len(ns)] for i in range(k)]
    if any(nn <= 0 for nn in ns) or any(xx < 0 for xx in xs) or \
            any(xx > nn for xx, nn in zip(xs, ns)):
        raise ValueError("Need n > 0 and 0 <= x <= n.")
    conf = _check_conf(conf)
    z = _qnorm_upper(conf)
    z2 = z * z  # R computes z^2 as z * z
    out = []
    for xx, nn in zip(xs, ns):
        centre = (xx + z2 / 2) / (nn + z2)
        half = z * math.sqrt(xx * (nn - xx) / nn + z2 / 4) / (nn + z2)
        out.append({"est": xx / nn, "lower": centre - half,
                    "upper": centre + half})
    return out[0] if scalar else out


def _quantile7(values: Sequence[float], probs: Sequence[float]) -> List[float]:
    """R's ``quantile(x, probs, type = 7)`` for finite ``values``."""
    x = sorted(values)
    n = len(x)
    out = []
    for pr in probs:
        index = 1 + max(n - 1, 0) * pr        # 1-based, as in R
        lo = math.floor(index)
        hi = math.ceil(index)
        qs = x[lo - 1]
        if index > lo and x[hi - 1] != qs:
            h = index - lo
            qs = (1 - h) * qs + h * x[hi - 1]
        out.append(qs)
    return out


def _is_na(v: Any) -> bool:
    """None, NaN (float or numpy) and pandas.NA / NaT count as missing."""
    if v is None:
        return True
    if isinstance(v, (str, bytes, bytearray)):
        return False
    if isinstance(v, numbers.Real):
        return math.isnan(float(v))
    return str(v) in ("<NA>", "NaT")


def cluster_boot_prop(x: Sequence[Any], cluster: Sequence[Hashable],
                      B: int = 2000, seed: Optional[int] = None,
                      conf: float = 0.95) -> Dict[str, Any]:
    """Cluster bootstrap confidence interval for a proportion.

    Percentile bootstrap interval for the proportion of ``True`` values in
    ``x``, resampling whole clusters (for example, studies that report
    several CFA models) with replacement. Quantiles use R's default
    (type 7) definition.

    Parameters
    ----------
    x
        Booleans (or 0/1 values). ``None``/NaN values are dropped.
    cluster
        Cluster identifier, same length as ``x``. Pairs whose cluster is
        missing are dropped too.
    B
        Number of bootstrap samples.
    seed
        Optional seed for a private ``random.Random(seed)`` generator. The
        global ``random`` state is never touched.
    conf
        Confidence level. Default 0.95.

    Returns
    -------
    dict
        Keys ``est`` (mean of ``x``), ``lower``, ``upper``, ``B``, ``n``
        (number of non-missing values) and ``n_clusters``.

    Notes
    -----
    Python's random number generator differs from R's, so for the same
    ``seed`` the bootstrap interval is *not* bit-identical to the R
    package's; it agrees only up to Monte Carlo error. ``est``, ``n`` and
    ``n_clusters`` are identical.

    Examples
    --------
    >>> x = [True, True, False, True, False, False, True, True]
    >>> study = [1, 1, 2, 3, 3, 3, 4, 5]
    >>> res = cluster_boot_prop(x, study, B=500, seed=1)
    >>> res["est"], res["n"], res["n_clusters"]
    (0.625, 8, 5)
    """
    xs = list(_as_list(x))
    cl = list(_as_list(cluster))
    if len(xs) != len(cl):
        raise ValueError("`x` and `cluster` must have the same length.")
    _check_conf(conf)
    if isinstance(B, bool) or not isinstance(B, numbers.Integral) or B < 1:
        raise ValueError("`B` must be a positive whole number.")
    B = int(B)
    keep = [(float(v), c) for v, c in zip(xs, cl)
            if not _is_na(v) and not _is_na(c)]
    if not keep:
        raise ValueError("No non-missing values in `x`.")

    # Per-cluster totals; clusters ordered like R's factor levels (sorted)
    # when the identifiers are sortable, else by first appearance.
    succ: Dict[Hashable, float] = {}
    size: Dict[Hashable, int] = {}
    for v, c in keep:
        succ[c] = succ.get(c, 0.0) + v
        size[c] = size.get(c, 0) + 1
    try:
        levels = sorted(succ)
    except TypeError:
        levels = list(succ)
    s = [succ[c] for c in levels]
    z = [size[c] for c in levels]
    k = len(levels)

    rng = random.Random(seed)
    boot = []
    for _ in range(B):
        pick = [rng.randrange(k) for _ in range(k)]
        boot.append(math.fsum(s[i] for i in pick) / sum(z[i] for i in pick))
    alpha = (1 - conf) / 2
    lower, upper = _quantile7(boot, [alpha, 1 - alpha])
    vals = [v for v, _ in keep]
    return {"est": math.fsum(vals) / len(vals), "lower": lower,
            "upper": upper, "B": B, "n": len(vals), "n_clusters": k}
