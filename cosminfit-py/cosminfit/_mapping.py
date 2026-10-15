"""Degrees of freedom, number of moments and the SRMR-RMSEA mapping.

Port of ``R/df.R`` and ``R/mapping.R`` of the R package cosminfit.

Vector arguments follow R's arithmetic recycling: the result has the length
of the longest argument and a warning is given when a longer length is not
a multiple of a shorter one. When every argument is a scalar the result is
a scalar; otherwise it is a list. R's ``NA`` in a numeric result is
returned as ``float("nan")``.
"""

from __future__ import annotations

import math
import numbers
import warnings
from typing import Any, Dict, List, Union

from ._rules import _as_list, _is_vector, _recycle_frame

__all__ = ["cfa_df", "n_moments", "srmr_to_rmsea", "rmsea_to_srmr",
           "rmsea_equivalent_table"]

Number = Union[int, float]
NumOrList = Union[Number, List[Number]]
CONVENTIONS = ("lavaan", "offdiag")


def _check_convention(convention: Any) -> str:
    if not (isinstance(convention, str) and convention in CONVENTIONS):
        raise ValueError("`convention` must be one of: " + ", ".join(CONVENTIONS))
    return convention


def _check_whole(x: Any, name: str) -> List[int]:
    """Return ``x`` as a list of ints; error unless all positive whole numbers."""
    vals = _as_list(x)
    ok = len(vals) > 0
    out: List[int] = []
    for v in vals:
        if (isinstance(v, bool) or not isinstance(v, numbers.Real)
                or not math.isfinite(float(v)) or v < 1 or v != round(v)):
            ok = False
            break
        out.append(int(v))
    if not ok:
        raise ValueError(f"`{name}` must contain positive whole numbers.")
    return out


def _recycle_arith(*vecs: List[Any]) -> List[List[Any]]:
    """R arithmetic recycling, with R's warning for non-multiple lengths."""
    lens = [len(v) for v in vecs]
    if min(lens) == 0:
        return [[] for _ in vecs]
    n = max(lens)
    if any(n % ln for ln in lens):
        warnings.warn("longer object length is not a multiple of shorter "
                      "object length", RuntimeWarning, stacklevel=3)
    return [v if len(v) == n else [v[i % len(v)] for i in range(n)]
            for v in vecs]


def _as_float(v: Any, name: str) -> float:
    """Numeric input for the mapping functions; None -> NaN."""
    if v is None:
        return math.nan
    if isinstance(v, (str, bytes, bytearray)):
        raise TypeError(f"`{name}` must be numeric.")
    try:
        return float(v)
    except (TypeError, ValueError):
        # pandas.NA and similar missing markers
        if str(v) in ("<NA>", "NaT"):
            return math.nan
        raise TypeError(f"`{name}` must be numeric.") from None


def _out(vals: List[Any], scalar: bool) -> Any:
    return vals[0] if scalar and len(vals) == 1 else vals


def cfa_df(p: Any, m: Any = 1) -> Union[int, List[int]]:
    """Degrees of freedom of a simple-structure CFA.

    Degrees of freedom of a CFA model with ``p`` items loading on ``m``
    correlated factors, each item on exactly one factor:
    ``p * (p + 1) / 2 - (2 * p + m * (m - 1) / 2)``.

    Parameters
    ----------
    p
        Number of items (scalar or sequence).
    m
        Number of factors (scalar or sequence, recycled with ``p``).

    Returns
    -------
    int or list of int
        Degrees of freedom. A ``ValueError`` is raised if any df is
        negative. A df of 0 gives a ``UserWarning``, because the SRMR-RMSEA
        mapping is then undefined.

    Examples
    --------
    >>> cfa_df(4)
    2
    >>> cfa_df(10, m=2)
    34
    >>> cfa_df([4, 6, 8, 10])
    [2, 9, 20, 35]
    """
    scalar = not (_is_vector(p) or _is_vector(m))
    pv = _check_whole(p, "p")
    mv = _check_whole(m, "m")
    pv, mv = _recycle_arith(pv, mv)
    # exact integer arithmetic: p(p+1) and m(m-1) are always even
    df = [pp * (pp + 1) // 2 - (2 * pp + mm * (mm - 1) // 2)
          for pp, mm in zip(pv, mv)]
    if any(d < 0 for d in df):
        raise ValueError("The model is not identified (negative df) for some "
                         "p and m.")
    if any(d == 0 for d in df):
        warnings.warn("df = 0 for some p and m; the SRMR-RMSEA mapping is "
                      "undefined.", UserWarning, stacklevel=2)
    return _out(df, scalar)


def n_moments(p: Any, convention: str = "lavaan") -> Union[int, List[int]]:
    """Number of sample moments ``D`` averaged by SRMR for ``p`` items.

    * ``"lavaan"``: ``D = p (p + 1) / 2``, including the diagonal. This is
      the convention of lavaan and of Shi, Maydeu-Olivares and Rosseel
      (2020), who divide by p(p+1)/2 even for correlation structures.
    * ``"offdiag"``: ``D = p (p - 1) / 2``, correlations below the diagonal
      only (an alternative, for sensitivity checks).

    Examples
    --------
    >>> n_moments(10)
    55
    >>> n_moments(10, convention="offdiag")
    45
    """
    convention = _check_convention(convention)
    pv = _check_whole(p, "p")
    if convention == "lavaan":
        d = [pp * (pp + 1) // 2 for pp in pv]
    else:
        d = [pp * (pp - 1) // 2 for pp in pv]
    return _out(d, not _is_vector(p))


def _positive_df(df: List[float]) -> List[float]:
    """Replace non-positive df by NaN, with a warning."""
    bad = [not math.isnan(d) and d <= 0 for d in df]
    if any(bad):
        warnings.warn("df must be positive; returning NA where df <= 0.",
                      UserWarning, stacklevel=4)
        df = [math.nan if b else d for d, b in zip(df, bad)]
    return df


def _convert(value: Any, name: str, p: Any, df: Any, convention: str,
             to_rmsea: bool) -> Union[float, List[float]]:
    convention = _check_convention(convention)
    scalar = not any(_is_vector(v) for v in (value, p, df))
    vals = [_as_float(v, name) for v in _as_list(value)]
    d = _as_list(n_moments(p, convention))
    dfs = _positive_df([_as_float(v, "df") for v in _as_list(df)])
    d, dfs = _recycle_arith(d, dfs)
    if to_rmsea:   # srmr * sqrt(D / df), operation order as in R
        ratio = [math.sqrt(a / b) if not math.isnan(b) else math.nan
                 for a, b in zip(d, dfs)]
    else:          # rmsea * sqrt(df / D)
        ratio = [math.sqrt(b / a) if not math.isnan(b) else math.nan
                 for a, b in zip(d, dfs)]
    vals, ratio = _recycle_arith(vals, ratio)
    res = [v * r for v, r in zip(vals, ratio)]
    return _out(res, scalar)


def srmr_to_rmsea(srmr: Any, p: Any, df: Any,
                  convention: str = "lavaan") -> Union[float, List[float]]:
    """Population RMSEA that corresponds to an SRMR value.

    Let ``G`` be the sum of squared residual correlations. Then
    ``SRMR = sqrt(G / D)`` and, at population level under unweighted least
    squares (ULS), ``RMSEA = sqrt(G / df)``. Hence
    ``RMSEA = SRMR * sqrt(D / df)`` and ``SRMR = RMSEA * sqrt(df / D)``,
    where ``D`` is given by :func:`n_moments`.

    This holds for population values under ULS. It does not hold exactly for
    sample values or for other estimators (such as ML), so do not use it to
    convert one reported index into another.

    Parameters
    ----------
    srmr
        SRMR value(s).
    p
        Number of items.
    df
        Model degrees of freedom (must be positive), e.g. from :func:`cfa_df`.
    convention
        Passed to :func:`n_moments`: ``"lavaan"`` (default) or ``"offdiag"``.

    Returns
    -------
    float or list of float
        NaN (with a warning) where ``df <= 0``.

    Examples
    --------
    >>> round(srmr_to_rmsea(0.08, p=4, df=2), 5)
    0.17889
    """
    return _convert(srmr, "srmr", p, df, convention, to_rmsea=True)


def rmsea_to_srmr(rmsea: Any, p: Any, df: Any,
                  convention: str = "lavaan") -> Union[float, List[float]]:
    """Population SRMR that corresponds to an RMSEA value.

    The inverse of :func:`srmr_to_rmsea`: ``SRMR = RMSEA * sqrt(df / D)``.

    Examples
    --------
    >>> round(rmsea_to_srmr(0.06, p=10, df=cfa_df(10)), 5)
    0.04786
    """
    return _convert(rmsea, "rmsea", p, df, convention, to_rmsea=False)


def rmsea_equivalent_table(p: Any = (4, 5, 6, 8, 10, 15, 20, 30), m: Any = 1,
                           srmr_cut: Any = 0.08,
                           convention: str = "lavaan") -> List[Dict[str, Any]]:
    """RMSEA equivalent of an SRMR cut-off, by number of items.

    For each number of items, the RMSEA value that corresponds, at
    population level under ULS, to the SRMR cut-off in a one- or
    multi-factor simple-structure CFA. See :func:`srmr_to_rmsea`.

    Parameters
    ----------
    p
        Numbers of items.
    m
        Number of factors (recycled with ``p``).
    srmr_cut
        SRMR cut-off. Default 0.08 (COSMIN).
    convention
        Passed to :func:`n_moments`.

    Returns
    -------
    list of dict
        One dict per row with keys ``p``, ``m``, ``df``, ``D`` and
        ``rmsea_equivalent`` (the R data.frame, row by row).

    Examples
    --------
    >>> [round(r["rmsea_equivalent"], 4) for r in rmsea_equivalent_table(p=[4, 30])]
    [0.1789, 0.0857]
    """
    convention = _check_convention(convention)
    n, cols = _recycle_frame({"p": _as_list(p), "m": _as_list(m)})
    pv, mv = cols["p"], cols["m"]
    df = _as_list(cfa_df(pv, mv))
    d = _as_list(n_moments(pv, convention))
    eq = _as_list(srmr_to_rmsea(srmr_cut, pv, df, convention))
    if len(eq) != n:
        raise ValueError(f"replacement has {len(eq)} rows, data has {n}")
    return [{"p": a, "m": b, "df": c, "D": e, "rmsea_equivalent": f}
            for a, b, c, e, f in zip(pv, mv, df, d, eq)]
