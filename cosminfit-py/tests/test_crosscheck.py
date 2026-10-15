"""Cross-language agreement with the R package.

Reads the CSV files written by ``crosscheck/make_vectors.R`` and requires
exact equality of every rating and |diff| < 1e-12 for every number.

The directory is ``$COSMINFIT_CROSSCHECK_DIR`` if set, else
``<repo>/crosscheck`` (two levels above this file). The tests are skipped
when the CSV files are not there.
"""

from __future__ import annotations

import csv
import math
import os
import sys
import unittest
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import cosminfit as cf
from cosminfit import _ci

HERE = Path(__file__).resolve().parent
XDIR = Path(os.environ.get("COSMINFIT_CROSSCHECK_DIR",
                           HERE.parent.parent / "crosscheck"))
TOL = 1e-12

#: agreement counts collected by the tests, printed at the end
SUMMARY: "dict[str, list[int]]" = defaultdict(lambda: [0, 0, 0])  # agree, total, bit-exact


def read_csv(name: str):
    path = XDIR / name
    if not path.exists():
        raise unittest.SkipTest(f"{path} not found; run crosscheck/make_vectors.R")
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def num(s: str):
    """CSV cell -> float or None (R's NA)."""
    return None if s == "NA" else float(s)


def numlist(s: str):
    return [num(v) for v in s.split(";")]


def record(key: str, ok: bool, exact: bool = True) -> None:
    c = SUMMARY[key]
    c[0] += ok
    c[1] += 1
    c[2] += exact


def close(a: float, b: float) -> "tuple[bool, bool]":
    if math.isnan(a) and math.isnan(b):
        return True, True
    return abs(a - b) < TOL, a == b


class TestRateCfaVectors(unittest.TestCase):
    """test_vectors.csv: rate_cfa() ratings, 8 rules x 2 missing modes."""

    @classmethod
    def setUpClass(cls):
        rows = read_csv("test_vectors.csv")
        cls.groups = defaultdict(list)
        for r in rows:
            cls.groups[(r["set"], r["rule"], r["missing"])].append(r)
        cls.n_rows = len(rows)

    def _check(self, as_strings: bool):
        mismatches = []
        for (set_, rule, mode), rows in sorted(self.groups.items()):
            if not as_strings and set_ == "d":
                continue  # set d *is* the string-parsing test
            cols = {}
            for k in ("cfi", "tli", "rmsea", "srmr", "p", "df"):
                vals = [r[k] for r in rows]
                if not as_strings or k in ("p", "df"):
                    vals = [num(v) for v in vals]
                cols[k] = vals
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                got = cf.rate_cfa(rule=rule, missing=mode, **cols)
            label = f"rate_cfa set {set_} ({'string' if as_strings else 'float'} input)"
            for r, g in zip(rows, got):
                ok = g == r["rating"]
                record(label, ok)
                if not ok:
                    mismatches.append((set_, r["case"], rule, mode, r["rating"], g))
        self.assertEqual(mismatches, [], f"{len(mismatches)} mismatches, first: "
                                         f"{mismatches[:5]}")

    def test_string_input(self):
        """CSV cells passed as strings (exercises the input cleaning)."""
        self._check(as_strings=True)

    def test_float_input(self):
        """CSV cells converted to float / None first."""
        self._check(as_strings=False)

    def test_vector_counts(self):
        sets = Counter(k[0] for k, v in self.groups.items() for _ in v)
        self.assertEqual(sets["a"], 81 * 16)
        self.assertEqual(sets["b"], 3000 * 16)
        self.assertGreater(sets["c"], 0)
        self.assertGreater(sets["d"], 0)
        self.assertEqual(len({k[1] for k in self.groups}), 8)
        self.assertEqual({k[2] for k in self.groups}, {"reported", "unknown"})


class TestRerateVectors(unittest.TestCase):
    """rerate_vectors.csv: rerate() with all 8 rules, n_indices, rmsea_rescue."""

    @classmethod
    def setUpClass(cls):
        cls.rows = read_csv("rerate_vectors.csv")
        cls.rules = [r["rule"] for r in cf.cosmin_rules()]

    def _compare(self, label, rows, out):
        bad = 0
        for r, o in zip(rows, out):
            want_resc = {"TRUE": True, "FALSE": False, "NA": None}[r["rmsea_rescue"]]
            got_resc = o["rmsea_rescue"]
            if got_resc is not None and not isinstance(got_resc, bool):
                got_resc = None if str(got_resc) == "<NA>" else bool(got_resc)
            ok = (all(o["rating_" + ru] == r["rating_" + ru] for ru in self.rules)
                  and int(o["n_indices"]) == int(r["n_indices"])
                  and got_resc == want_resc)
            record(label, ok)
            bad += not ok
        self.assertEqual(bad, 0)
        self.assertEqual(len(out), len(rows))

    def test_list_of_dicts(self):
        for mode in ("reported", "unknown"):
            rows = [r for r in self.rows if r["missing"] == mode]
            data = [{k: r[k] for k in ("cfi", "tli", "rmsea", "srmr", "p", "df")}
                    for r in rows]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out = cf.rerate(data, rules=self.rules, missing=mode)
            self._compare("rerate (list of dicts)", rows, out)

    def test_pandas(self):
        try:
            import pandas as pd
        except ImportError:
            self.skipTest("pandas not installed")
        for mode in ("reported", "unknown"):
            rows = [r for r in self.rows if r["missing"] == mode]
            df = pd.DataFrame({k: [num(r[k]) for r in rows]
                               for k in ("cfi", "tli", "rmsea", "srmr", "p", "df")},
                              dtype=float)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                out = cf.rerate(df, rules=self.rules, missing=mode)
            recs = [{c: out[c].iloc[i] for c in out.columns} for i in range(len(out))]
            self._compare("rerate (pandas DataFrame)", rows, recs)


class TestMappingVectors(unittest.TestCase):
    """mapping_vectors.csv: df, moments, SRMR-RMSEA mapping, CIs."""

    @classmethod
    def setUpClass(cls):
        cls.by_func = defaultdict(list)
        for r in read_csv("mapping_vectors.csv"):
            cls.by_func[r["func"]].append(r)

    def _num_check(self, label, got, want):
        ok, exact = close(float(got), want)
        record(label, ok, exact)
        self.assertTrue(ok, f"{label}: got {got!r}, R {want!r}")

    def test_cfa_df(self):
        for r in self.by_func["cfa_df"]:
            self._num_check("cfa_df", cf.cfa_df(int(num(r["p"])), int(num(r["m"]))),
                            num(r["value"]))

    def test_n_moments(self):
        for r in self.by_func["n_moments"]:
            self._num_check("n_moments", cf.n_moments(int(num(r["p"])), r["convention"]),
                            num(r["value"]))

    def test_srmr_rmsea(self):
        for func, f in (("srmr_to_rmsea", cf.srmr_to_rmsea),
                        ("rmsea_to_srmr", cf.rmsea_to_srmr)):
            rows = self.by_func[func]
            self.assertGreater(len(rows), 0)
            for r in rows:   # scalar calls
                got = f(num(r["x"]), int(num(r["p"])), num(r["df"]), r["convention"])
                self._num_check(func, got, num(r["value"]))
            # one vectorised call per (x, convention)
            groups = defaultdict(list)
            for r in rows:
                groups[(r["x"], r["convention"])].append(r)
            for (x, conv), g in groups.items():
                got = f(num(x), [int(num(r["p"])) for r in g],
                        [num(r["df"]) for r in g], conv)
                for gv, r in zip(got, g):
                    self._num_check(func + " (vectorised)", gv, num(r["value"]))

    def test_rmsea_equivalent_table(self):
        groups = defaultdict(list)
        for r in self.by_func["rmsea_equivalent_table"]:
            groups[(r["convention"], r["m"], r["output"])].append(r)
        self.assertGreater(len(groups), 0)
        for (conv, m, output), g in groups.items():
            tab = cf.rmsea_equivalent_table(p=[int(num(r["p"])) for r in g],
                                            m=int(num(m)), convention=conv)
            self.assertEqual(len(tab), len(g))
            for row, r in zip(tab, g):
                self.assertEqual(row["df"], num(r["df"]))
                self._num_check("rmsea_equivalent_table", row[output], num(r["value"]))

    def test_wilson_ci(self):
        for r in self.by_func["wilson_ci"]:
            w = cf.wilson_ci(num(r["x"]), num(r["n"]), conf=num(r["conf"]))
            self._num_check("wilson_ci", w[r["output"]], num(r["value"]))

    def test_qnorm(self):
        for r in self.by_func["qnorm"]:
            self._num_check("qnorm (statistics.NormalDist)",
                            _ci._qnorm_upper(num(r["conf"])), num(r["value"]))

    def test_quantile_type7(self):
        for r in self.by_func["quantile_type7"]:
            got = _ci._quantile7(numlist(r["x"]), [num(r["conf"])])[0]
            self._num_check("quantile type 7", got, num(r["value"]))

    def test_cluster_boot_prop(self):
        rows = {r["output"]: r for r in self.by_func["cluster_boot_prop"]}
        self.assertTrue(rows)
        r0 = rows["est"]
        x = [None if v is None else bool(v) for v in numlist(r0["x"])]
        cl = [int(v) for v in numlist(r0["n"])]
        res = cf.cluster_boot_prop(x, cl, B=20000, seed=7)
        # deterministic parts: exact
        for o in ("est", "B", "n", "n_clusters"):
            self._num_check("cluster_boot_prop (deterministic parts)", res[o],
                            num(rows[o]["value"]))
        # interval: different RNG, so only Monte Carlo agreement
        for o in ("lower", "upper"):
            diff = abs(res[o] - num(rows[o]["value"]))
            record("cluster_boot_prop CI within 0.02 of R (MC only)", diff < 0.02,
                   exact=False)
            self.assertLess(diff, 0.02)


def tearDownModule():
    if not SUMMARY:
        return
    lines = ["", "Cross-check against R (agree / total; bit-identical numbers):"]
    for k in sorted(SUMMARY):
        a, t, e = SUMMARY[k]
        pct = 100.0 * a / t if t else 0.0
        lines.append(f"  {k:<52s} {a:>6d} / {t:<6d} ({pct:6.2f}%)  exact {e}")
    tot_a = sum(v[0] for v in SUMMARY.values())
    tot_t = sum(v[1] for v in SUMMARY.values())
    lines.append(f"  {'TOTAL':<52s} {tot_a:>6d} / {tot_t:<6d} "
                 f"({100.0 * tot_a / tot_t:6.2f}%)")
    sys.stderr.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    unittest.main()
