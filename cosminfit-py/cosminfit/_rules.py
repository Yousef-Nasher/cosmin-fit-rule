"""Rating rules, input cleaning and the rule-tree engine.

This is a line-by-line port of ``R/rules.R`` of the R package cosminfit.
Every rule is a small logic tree over "pass indicators" (leaves such as
``"cfi_gt95"``). Nodes are ``("or", args)``, ``("and", args)`` and
``("at_least", k, args)``; leaves are indicator names.
"""

from __future__ import annotations

import itertools
import math
import re
import warnings
from collections.abc import Iterable, Mapping
from typing import Any, Dict, List, Optional, Tuple, Union

__all__ = ["cosmin_rules", "rate_cfa"]

Rating = str
Tri = Optional[bool]  # True / False / None (= not reported / not evaluable)
Node = Union[str, Tuple[Any, ...]]

MISSING_MODES = ("reported", "unknown")


def _or(*args: Node) -> Node:
    return ("or", args)


def _and(*args: Node) -> Node:
    return ("and", args)


def _at_least(k: int, *args: Node) -> Node:
    return ("at_least", k, args)


RULE_TREES: Dict[str, Node] = {
    "cosmin": _or("cfi_gt95", "tli_gt95", "rmsea_lt06", "srmr_lt08"),
    "no_rmsea": _or("cfi_gt95", "tli_gt95", "srmr_lt08"),
    "srmr_df": _or("cfi_gt95", "tli_gt95", "rmsea_lt06", "srmr_lt_df"),
    "two_of_three": _at_least(2, "cfi_gt95", "rmsea_lt06", "srmr_lt08"),
    "all_three": _and("cfi_gt95", "rmsea_lt06", "srmr_lt08"),
    "hu_bentler": _or(_and("cfi_ge95", "srmr_le09"),
                      _and("rmsea_le06", "srmr_le09")),
    "prinsen2016_a": _and(_or("cfi_gt95", "tli_gt95"),
                          _or("rmsea_lt06", "srmr_lt08")),
    "prinsen2016_b": _or(_and(_or("cfi_gt95", "tli_gt95"), "rmsea_lt06"),
                         "srmr_lt08"),
}

_DEFINITIONS = [
    "CFI > .95 OR TLI > .95 OR RMSEA < .06 OR SRMR < .08",
    "CFI > .95 OR TLI > .95 OR SRMR < .08",
    "CFI > .95 OR TLI > .95 OR RMSEA < .06 OR SRMR < .06*sqrt(df/D)",
    "at least two of: CFI > .95, RMSEA < .06, SRMR < .08",
    "CFI > .95 AND RMSEA < .06 AND SRMR < .08",
    "(CFI >= .95 AND SRMR <= .09) OR (RMSEA <= .06 AND SRMR <= .09)",
    "(CFI > .95 OR TLI > .95) AND (RMSEA < .06 OR SRMR < .08)",
    "((CFI > .95 OR TLI > .95) AND RMSEA < .06) OR SRMR < .08",
]
_INDICES = ["cfi, tli, rmsea, srmr", "cfi, tli, srmr", "cfi, tli, rmsea, srmr",
            "cfi, rmsea, srmr", "cfi, rmsea, srmr", "cfi, rmsea, srmr",
            "cfi, tli, rmsea, srmr", "cfi, tli, rmsea, srmr"]
_NEEDS_P_DF = [False, False, True, False, False, False, False, False]


def cosmin_rules() -> List[Dict[str, Any]]:
    """List the rules that :func:`rate_cfa` and :func:`rerate` can apply.

    All inequalities are as written (strict unless shown otherwise).

    * ``cosmin``: CFI > .95 OR TLI > .95 OR RMSEA < .06 OR SRMR < .08
      (Prinsen et al. 2018; Mokkink et al. 2024).
    * ``no_rmsea``: CFI > .95 OR TLI > .95 OR SRMR < .08.
    * ``srmr_df``: CFI > .95 OR TLI > .95 OR RMSEA < .06 OR
      SRMR < .06 * sqrt(df / D), with D = p(p+1)/2. Needs ``p`` and ``df``.
    * ``two_of_three``: at least two of CFI > .95, RMSEA < .06, SRMR < .08
      (TLI not used).
    * ``all_three``: CFI > .95 AND RMSEA < .06 AND SRMR < .08.
    * ``hu_bentler``: (CFI >= .95 AND SRMR <= .09) OR
      (RMSEA <= .06 AND SRMR <= .09).
    * ``prinsen2016_a``, ``prinsen2016_b``: two readings of the 2016
      COSMIN-COMET guideline (Prinsen et al. 2016, Trials 17:449, Table 4).
      Reading a is (CFI > .95 OR TLI > .95) AND (RMSEA < .06 OR SRMR < .08);
      reading b is ((CFI > .95 OR TLI > .95) AND RMSEA < .06) OR SRMR < .08.

    Returns
    -------
    list of dict
        One dict per rule with keys ``rule``, ``definition``, ``indices``
        (the indices the rule uses) and ``needs_p_df`` (whether ``p`` and
        ``df`` are needed). This mirrors the R data.frame row by row.
    """
    return [
        {"rule": r, "definition": d, "indices": i, "needs_p_df": n}
        for r, d, i, n in zip(RULE_TREES, _DEFINITIONS, _INDICES, _NEEDS_P_DF)
    ]


# ---- input helpers ----------------------------------------------------------

_R_SPACE = " \t\n\v\f\r"
_DEC_RE = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d*)?")
_DEC_EXP_EMPTY = re.compile(r"[eE][+-]?$")
_HEX_RE = re.compile(r"[+-]?0[xX](?:[0-9a-fA-F]+\.?[0-9a-fA-F]*|\.[0-9a-fA-F]+)"
                     r"(?:[pP][+-]?\d+)?")


def _parse_string(s: str) -> Optional[float]:
    """Read a string as a number the way R's ``as.numeric()`` does.

    Leading/trailing ASCII white space is ignored; decimal and hexadecimal
    (``0x1A``) forms are accepted, as is a dangling exponent (``"1e"``, which
    R reads as 1). Python-only forms (``"1_000"``, non-ASCII digits) and
    anything else give ``None``.
    """
    t = s.strip(_R_SPACE)
    if not t or not t.isascii():
        return None
    if _DEC_RE.fullmatch(t):
        t = _DEC_EXP_EMPTY.sub("", t)
        return float(t)
    if _HEX_RE.fullmatch(t):
        try:
            return float.fromhex(t)
        except (ValueError, OverflowError):
            return None
    return None  # "NA", "NaN", "Inf", "n.r." ... are all "not reported"


def _as_index_value(v: Any) -> Optional[float]:
    """Coerce one value to a finite float, or ``None`` (= not reported)."""
    if v is None:
        return None
    if isinstance(v, str):
        x = _parse_string(v)
    elif isinstance(v, (bytes, bytearray)):
        try:
            x = _parse_string(bytes(v).decode("ascii"))
        except UnicodeDecodeError:
            x = None
    else:
        try:
            x = float(v)
        except (TypeError, ValueError, OverflowError):
            # e.g. pandas.NA, pandas.NaT, objects without __float__
            x = None
    if x is None or not math.isfinite(x):
        return None
    return x


def _is_vector(x: Any) -> bool:
    """True for list-like input (list, tuple, range, Series, ndarray...)."""
    if isinstance(x, (str, bytes, bytearray)):
        return False
    if getattr(x, "ndim", None) == 0:  # numpy 0-d array / numpy scalar
        return False
    return isinstance(x, Iterable) and hasattr(x, "__len__")


def _as_list(x: Any) -> List[Any]:
    """Scalar -> one-element list; Mapping -> its values; else list(x)."""
    if not _is_vector(x):
        return [x]
    if isinstance(x, Mapping):  # like an R named vector: names are dropped
        return list(x.values())
    return list(x)


def as_index(x: Any) -> List[Optional[float]]:
    """Vector version of the R helper ``as_index()``.

    Anything that is not a finite number (None, NaN, +/-inf, empty or
    non-numeric strings, pandas.NA) becomes ``None``. Numeric strings are
    parsed.
    """
    return [_as_index_value(v) for v in _as_list(x)]


def _recycle_frame(cols: Dict[str, List[Any]]) -> Tuple[int, Dict[str, List[Any]]]:
    """Recycle columns like ``as.data.frame(list(...))`` in R.

    Every length must divide the longest one; otherwise a ``ValueError``.
    """
    lens = {k: len(v) for k, v in cols.items()}
    n = max(lens.values()) if lens else 0
    if n == 0:
        return 0, {k: [] for k in cols}
    bad = sorted({ln for ln in lens.values() if ln == 0 or n % ln != 0})
    if bad:
        shown = ", ".join(str(v) for v in sorted(set(lens.values())))
        raise ValueError(f"arguments imply differing number of rows: {shown}")
    out = {}
    for k, v in cols.items():
        out[k] = v if len(v) == n else [v[i % len(v)] for i in range(n)]
    return n, out


def check_rule(rule: Any) -> str:
    """Raise ``ValueError`` unless ``rule`` is one known rule name."""
    if not (isinstance(rule, str) and rule in RULE_TREES):
        raise ValueError("`rule` must be one of: " + ", ".join(RULE_TREES))
    return rule


def check_missing(missing: Any) -> str:
    """Raise ``ValueError`` unless ``missing`` is "reported" or "unknown"."""
    if not (isinstance(missing, str) and missing in MISSING_MODES):
        raise ValueError("`missing` must be one of: " + ", ".join(MISSING_MODES))
    return missing


def check_ranges(x: Dict[str, List[Optional[float]]]) -> None:
    """Warn about implausible index values (as in the R package)."""
    if any(v is not None and v > 1.2 for v in x["cfi"] + x["tli"]):
        warnings.warn("Some CFI or TLI values are above 1.2; check the input.",
                      UserWarning, stacklevel=3)
    if any(v is not None and v < 0 for v in x["rmsea"] + x["srmr"]):
        warnings.warn("Some RMSEA or SRMR values are negative; check the input.",
                      UserWarning, stacklevel=3)


def _cmp(a: Optional[float], op: str, b: Optional[float]) -> Tri:
    if a is None or b is None:
        return None
    if op == ">":
        return a > b
    if op == ">=":
        return a >= b
    if op == "<":
        return a < b
    return a <= b


def pass_indicators(x: Dict[str, List[Optional[float]]]) -> List[Dict[str, Tri]]:
    """Per-row pass indicators; ``None`` means "not reported"."""
    rows = []
    for cfi, tli, rmsea, srmr, p, df in zip(x["cfi"], x["tli"], x["rmsea"],
                                            x["srmr"], x["p"], x["df"]):
        cut: Optional[float] = None
        if p is not None and df is not None and p >= 1 and df > 0:
            # same operation order as R: 0.06 * sqrt(df / (p * (p + 1) / 2))
            cut = 0.06 * math.sqrt(df / (p * (p + 1) / 2))
        rows.append({
            "cfi_gt95": _cmp(cfi, ">", 0.95),
            "tli_gt95": _cmp(tli, ">", 0.95),
            "rmsea_lt06": _cmp(rmsea, "<", 0.06),
            "srmr_lt08": _cmp(srmr, "<", 0.08),
            "srmr_lt_df": _cmp(srmr, "<", cut),
            "cfi_ge95": _cmp(cfi, ">=", 0.95),
            "rmsea_le06": _cmp(rmsea, "<=", 0.06),
            "srmr_le09": _cmp(srmr, "<=", 0.09),
        })
    return rows


# ---- engine -----------------------------------------------------------------

def eval_reported(node: Node, ind: Mapping[str, Tri]) -> Tri:
    """Evaluate a tree with "reported" semantics; ``None`` = not evaluable."""
    if isinstance(node, str):
        return ind[node]
    op = node[0]
    args = node[2] if op == "at_least" else node[1]
    vals = [eval_reported(a, ind) for a in args]
    n_eval = sum(v is not None for v in vals)
    n_true = sum(v is True for v in vals)
    if op == "or":
        return (n_true > 0) if n_eval > 0 else None
    if op == "and":
        return (n_true == len(vals)) if n_eval == len(vals) else None
    k = node[1]
    return (n_true >= k) if n_eval >= k else None


def tree_leaves(node: Node) -> List[str]:
    """Leaf names of a tree, depth first (may contain duplicates)."""
    if isinstance(node, str):
        return [node]
    args = node[2] if node[0] == "at_least" else node[1]
    return [leaf for a in args for leaf in tree_leaves(a)]


def eval_unknown(node: Node, ind: Mapping[str, Tri]) -> Tri:
    """Evaluate with "unknown" semantics: try every completion of the
    missing leaves; keep the result only if it is the same for all."""
    leaves = list(dict.fromkeys(tree_leaves(node)))
    missing = [lf for lf in leaves if ind[lf] is None]
    results = set()
    for combo in itertools.product((False, True), repeat=len(missing)):
        filled = dict(ind)
        filled.update(zip(missing, combo))
        results.add(eval_reported(node, filled))
        if len(results) > 1:
            return None
    (res,) = results
    return res


def rate_core(ind: List[Dict[str, Tri]], rule: str, missing: str) -> List[Rating]:
    """Rate a list of indicator rows under one rule."""
    tree = RULE_TREES[rule]
    ev = eval_reported if missing == "reported" else eval_unknown
    out = []
    for row in ind:
        res = ev(tree, row)
        out.append("sufficient" if res is True else
                   "insufficient" if res is False else "indeterminate")
    return out


def rate_cfa(cfi: Any = None, tli: Any = None, rmsea: Any = None,
             srmr: Any = None, rule: str = "cosmin", p: Any = None,
             df: Any = None, missing: str = "reported"
             ) -> Union[Rating, List[Rating]]:
    """Rate CFA fit against COSMIN and alternative rules.

    Rates each set of fit indices as ``"sufficient"``, ``"insufficient"``
    or ``"indeterminate"`` under one rule (see :func:`cosmin_rules`).
    Inputs may be scalars or sequences; sequences are recycled like
    data.frame columns in R (every length must divide the longest).
    Values that are ``None``, NaN, infinite, or cannot be read as a number
    count as "not reported"; numeric strings such as ``"0.97"`` are parsed.

    Each rule is a logic formula over pass indicators (for example
    CFI > .95). An indicator is missing when its index is not reported. For
    ``srmr_df`` the SRMR term is missing also when ``p`` or ``df`` is
    missing, ``p < 1`` or ``df <= 0``.

    ``missing="reported"`` (primary analysis) evaluates the formula
    bottom-up: an OR group is evaluable if at least one member is evaluable
    and is True if any evaluable member is True; an AND group is evaluable
    only if all members are; "at least k of" is evaluable if at least k
    members are, and True if at least k evaluable members are True.
    Not evaluable gives ``"indeterminate"``.

    ``missing="unknown"`` (sensitivity analysis) lets each missing indicator
    be True or False, enumerates all completions, and gives a definite
    rating only if every completion agrees (three-valued Kleene logic).

    With all indices of a rule reported both modes give the same rating.

    Parameters
    ----------
    cfi, tli, rmsea, srmr
        Fit index values.
    rule
        One rule name from :func:`cosmin_rules`. Default ``"cosmin"``.
    p
        Number of items (only used by ``srmr_df``).
    df
        Model degrees of freedom (only used by ``srmr_df``).
    missing
        ``"reported"`` (default) or ``"unknown"``.

    Returns
    -------
    str or list of str
        A single rating when every argument is a scalar, otherwise a list
        (empty if any argument has length zero).

    Examples
    --------
    >>> rate_cfa(cfi=0.93, rmsea=0.05)
    'sufficient'
    >>> rate_cfa(cfi=0.97, rule="all_three")
    'indeterminate'
    >>> rate_cfa(cfi=[0.96, 0.90, None], srmr=[None, 0.09, None])
    ['sufficient', 'insufficient', 'indeterminate']
    """
    missing = check_missing(missing)
    rule = check_rule(rule)
    args = {"cfi": cfi, "tli": tli, "rmsea": rmsea, "srmr": srmr,
            "p": p, "df": df}
    scalar = not any(_is_vector(v) for v in args.values())
    cols = {k: as_index(v) for k, v in args.items()}
    if any(len(v) == 0 for v in cols.values()):
        return []
    _, x = _recycle_frame(cols)
    check_ranges(x)
    out = rate_core(pass_indicators(x), rule, missing)
    return out[0] if scalar else out
