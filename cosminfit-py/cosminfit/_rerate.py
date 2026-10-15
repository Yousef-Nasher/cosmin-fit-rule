"""Re-rate a data set of CFA results under several rules (port of ``R/rerate.R``)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Dict, List, Optional

from ._rules import (as_index, check_missing, check_ranges, check_rule,
                     pass_indicators, rate_core)

__all__ = ["rerate"]

DEFAULT_RULES = ("cosmin", "no_rmsea", "srmr_df", "two_of_three", "all_three")
DEFAULT_COLS = {"cfi": "cfi", "tli": "tli", "rmsea": "rmsea", "srmr": "srmr",
                "p": "p", "df": "df"}


def _is_dataframe(data: Any) -> bool:
    try:
        import pandas as pd  # optional dependency
    except ImportError:
        return False
    return isinstance(data, pd.DataFrame)


def rerate(data: Any, rules: Sequence[str] = DEFAULT_RULES,
           missing: str = "reported",
           cols: Optional[Mapping[str, Optional[str]]] = None) -> Any:
    """Re-rate a data set of CFA results under several rules.

    Applies :func:`rate_cfa` with several rules to each row of ``data``.

    Parameters
    ----------
    data
        A pandas DataFrame, or a list of dicts (one dict per CFA result).
    rules
        Rule names from :func:`cosmin_rules` (a single string is allowed).
    missing
        ``"reported"`` (default) or ``"unknown"``; see :func:`rate_cfa`.
    cols
        Mapping giving the column names for ``cfi``, ``tli``, ``rmsea``,
        ``srmr``, ``p`` and ``df``. Partial mappings are filled in with the
        defaults (the same names). A column that is absent from ``data``
        (or a key missing from a row dict) counts as not reported.

    Returns
    -------
    pandas.DataFrame or list of dict
        A copy of ``data`` (same type) with one column ``rating_<rule>`` per
        rule, plus ``n_indices`` (how many of CFI, TLI, RMSEA and SRMR were
        reported) and ``rmsea_rescue`` (True when RMSEA < .06 is the only
        one of the four COSMIN clauses that passes among reported indices;
        ``None`` (``pandas.NA`` in a DataFrame) when RMSEA was not
        reported). The input is not modified.

    Examples
    --------
    >>> d = [{"cfi": 0.97, "tli": None, "rmsea": 0.05, "srmr": 0.04, "p": 10, "df": 35},
    ...      {"cfi": 0.91, "tli": 0.89, "rmsea": 0.05, "srmr": 0.09, "p": 6, "df": 9}]
    >>> [r["rating_all_three"] for r in rerate(d)]
    ['sufficient', 'insufficient']
    >>> [r["rmsea_rescue"] for r in rerate(d)]
    [False, True]
    """
    is_df = _is_dataframe(data)
    if not is_df:
        ok = (isinstance(data, Sequence) and not isinstance(data, (str, bytes))
              and all(isinstance(r, Mapping) for r in data))
        if not ok:
            raise TypeError("`data` must be a pandas DataFrame or a list of "
                            "dicts (the R version requires a data.frame).")
    missing = check_missing(missing)
    if isinstance(rules, str):
        rules = [rules]
    rules = list(rules)
    for r in rules:
        check_rule(r)
    colmap: Dict[str, Optional[str]] = dict(DEFAULT_COLS)
    if cols is not None:
        colmap.update(dict(cols))

    nrow = len(data)

    def get_col(name: str) -> List[Optional[float]]:
        col = colmap.get(name)
        if is_df:
            if col is not None and col in data.columns:
                return as_index(data[col].tolist())
            return [None] * nrow
        if col is None:
            return [None] * nrow
        return as_index([row.get(col) for row in data])

    x = {name: get_col(name) for name in DEFAULT_COLS}
    check_ranges(x)
    ind = pass_indicators(x)

    new_cols: Dict[str, List[Any]] = {}
    for r in rules:
        new_cols["rating_" + r] = rate_core(ind, r, missing)
    new_cols["n_indices"] = [
        sum(v is not None for v in vals)
        for vals in zip(x["cfi"], x["tli"], x["rmsea"], x["srmr"])
    ]
    rescue: List[Optional[bool]] = []
    for row in ind:
        if row["rmsea_lt06"] is None:
            rescue.append(None)
        else:
            other_pass = any(row[k] is True
                             for k in ("cfi_gt95", "tli_gt95", "srmr_lt08"))
            rescue.append(row["rmsea_lt06"] and not other_pass)
    new_cols["rmsea_rescue"] = rescue

    if is_df:
        import pandas as pd
        out = data.copy()
        for k, v in new_cols.items():
            if k == "rmsea_rescue":
                out[k] = pd.array(v, dtype="boolean")
            else:
                out[k] = v
        return out
    result = []
    for i, row in enumerate(data):
        new = dict(row)
        for k, v in new_cols.items():
            new[k] = v[i]
        result.append(new)
    return result
