"""Port of the R package's unit tests (cosminfit/tests/unit/*.R).

Each R ``check()`` becomes an assertion; the order follows the R files.
"""

from __future__ import annotations

import itertools
import math
import random
import unittest
import warnings
from statistics import NormalDist

import cosminfit as cf

RULES = ["cosmin", "no_rmsea", "srmr_df", "two_of_three", "all_three",
         "hu_bentler", "prinsen2016_a", "prinsen2016_b"]


class WarnMixin:
    def assertWarns_(self, fn, *a, **k):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = fn(*a, **k)
        self.assertTrue(w, "expected a warning")
        return res

    def assertNoWarning(self, fn, *a, **k):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            res = fn(*a, **k)
        self.assertEqual([str(x.message) for x in w], [])
        return res

    def quiet(self, fn, *a, **k):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return fn(*a, **k)


# ---- test-ci.R ----------------------------------------------------------------

class TestCI(unittest.TestCase):
    def test_wilson(self):
        w = cf.wilson_ci(5, 100)
        z = NormalDist().inv_cdf(0.975)
        centre = (5 + z ** 2 / 2) / (100 + z ** 2)
        half = z / (100 + z ** 2) * math.sqrt(5 * 95 / 100 + z ** 2 / 4)
        self.assertLess(abs(w["est"] - 0.05), 1e-12)
        self.assertLess(abs(w["lower"] - (centre - half)), 1e-12)
        self.assertLess(abs(w["upper"] - (centre + half)), 1e-12)
        self.assertLess(abs(w["lower"] - 0.0215), 1e-4)   # hand values
        self.assertLess(abs(w["upper"] - 0.1118), 1e-4)
        # exact R values for wilson_ci(5, 100)
        self.assertEqual(w["lower"], 0.02154367915436798)   # R, sprintf("%.17g")
        self.assertEqual(w["upper"], 0.11175046923191911)
        w0 = cf.wilson_ci([0, 50], 50)
        self.assertEqual(len(w0), 2)
        self.assertLess(abs(w0[0]["lower"]), 1e-12)
        self.assertLess(abs(w0[1]["upper"] - 1), 1e-12)
        self.assertLess(cf.wilson_ci(5, 100, conf=0.99)["lower"], w["lower"])
        with self.assertRaises(ValueError):
            cf.wilson_ci(6, 5)
        with self.assertRaises(ValueError):
            cf.wilson_ci(1, 5, conf=1.2)
        with self.assertRaises(ValueError):
            cf.wilson_ci(0, 0)
        with self.assertRaises(ValueError):
            cf.wilson_ci(None, 5)

    def test_cluster_boot(self):
        rng = random.Random(42)
        x = [rng.random() < 0.4 for _ in range(60)]
        x[2] = None
        x[16] = None
        cl = [i // 3 + 1 for i in range(60)]
        b1 = cf.cluster_boot_prop(x, cl, B=500, seed=7)
        b2 = cf.cluster_boot_prop(x, cl, B=500, seed=7)
        self.assertEqual(b1, b2)                         # reproducible
        obs = [float(v) for v in x if v is not None]
        self.assertLess(abs(b1["est"] - sum(obs) / len(obs)), 1e-12)
        self.assertTrue(b1["lower"] <= b1["est"] <= b1["upper"])
        self.assertEqual((b1["n"], b1["n_clusters"], b1["B"]), (58, 20, 500))
        self.assertEqual(list(b1), ["est", "lower", "upper", "B", "n", "n_clusters"])
        # global RNG state untouched
        random.seed(99)
        r1 = random.random()
        random.seed(99)
        cf.cluster_boot_prop(x, cl, B=50, seed=1)
        r2 = random.random()
        self.assertEqual(r1, r2)
        with self.assertRaises(ValueError):
            cf.cluster_boot_prop(x, cl[1:])
        with self.assertRaises(ValueError):
            cf.cluster_boot_prop([None, None], [1, 2])
        with self.assertRaises(ValueError):
            cf.cluster_boot_prop(x, cl, conf=1)
        # NaN values and missing clusters are dropped too
        b3 = cf.cluster_boot_prop([1, 0, float("nan"), 1], [1, 2, 3, None], B=10, seed=1)
        self.assertEqual((b3["n"], b3["n_clusters"], b3["est"]), (2, 2, 0.5))


# ---- test-df-mapping.R ------------------------------------------------------------

class TestDfMapping(WarnMixin, unittest.TestCase):
    def test_cfa_df(self):
        self.assertEqual(cf.cfa_df(4, 1), 2)
        self.assertEqual(cf.cfa_df(10, 2), 34)
        self.assertEqual(cf.cfa_df(30, 1), 405)
        self.assertEqual(cf.cfa_df([4, 10, 30], [1, 2, 1]), [2, 34, 405])
        self.assertEqual(cf.cfa_df(4.0), 2)
        self.assertEqual(self.assertWarns_(cf.cfa_df, 3), 0)    # df = 0 warns
        for bad in (2, 4.5, 0, -1, True, "4", float("nan"), float("inf"), []):
            with self.assertRaises(ValueError, msg=repr(bad)):
                self.quiet(cf.cfa_df, bad)

    def test_n_moments(self):
        self.assertEqual(cf.n_moments(10), 55)
        self.assertEqual(cf.n_moments(10, "offdiag"), 45)
        self.assertEqual(cf.n_moments([4, 10]), [10, 55])
        with self.assertRaises(ValueError):
            cf.n_moments(10, "other")

    def test_rmsea_equivalent_table(self):
        tab = cf.rmsea_equivalent_table()
        self.assertEqual(list(tab[0]), ["p", "m", "df", "D", "rmsea_equivalent"])
        self.assertEqual(len(tab), 8)
        r4 = [r for r in tab if r["p"] == 4][0]
        self.assertEqual((r4["df"], r4["D"]), (2, 10))
        self.assertLess(abs(r4["rmsea_equivalent"] - 0.1789), 1e-3)
        self.assertEqual(round(r4["rmsea_equivalent"], 5), 0.17889)
        r30 = [r for r in tab if r["p"] == 30][0]
        self.assertEqual((r30["df"], r30["D"]), (405, 465))
        self.assertLess(abs(r30["rmsea_equivalent"] - 0.0857), 1e-3)
        self.assertEqual(round(r30["rmsea_equivalent"], 5), 0.08572)
        eq = [r["rmsea_equivalent"] for r in tab]
        self.assertTrue(all(b < a for a, b in zip(eq, eq[1:])))
        off = cf.rmsea_equivalent_table(p=4, convention="offdiag")
        self.assertLess(abs(off[0]["rmsea_equivalent"] - 0.08 * math.sqrt(6 / 2)), 1e-12)
        t3 = cf.rmsea_equivalent_table(p=[6, 9, 12], m=3)
        self.assertEqual([r["df"] for r in t3], [6, 24, 51])
        with self.assertRaises(ValueError):
            cf.rmsea_equivalent_table(p=[4, 5, 6], m=[1, 2])

    def test_srmr_rmsea(self):
        self.assertLess(abs(cf.srmr_to_rmsea(0.08, 4, 2) - 0.08 * math.sqrt(5)), 1e-12)
        s = [round(0.01 * i, 2) for i in range(1, 16)]
        for conv in ("lavaan", "offdiag"):
            back = cf.rmsea_to_srmr(cf.srmr_to_rmsea(s, 12, cf.cfa_df(12), conv),
                                    12, cf.cfa_df(12), conv)
            self.assertLess(max(abs(a - b) for a, b in zip(back, s)), 1e-12)
            back2 = cf.srmr_to_rmsea(cf.rmsea_to_srmr(s, 7, 14, conv), 7, 14, conv)
            self.assertLess(max(abs(a - b) for a, b in zip(back2, s)), 1e-12)
        v = self.assertWarns_(cf.srmr_to_rmsea, 0.05, 3, 0)    # df = 0 warns
        self.assertTrue(math.isnan(v))                           # ... gives NA
        self.assertTrue(math.isnan(cf.srmr_to_rmsea(None, 4, 2)))
        self.assertTrue(math.isnan(cf.srmr_to_rmsea(0.08, 4, None)))
        with self.assertRaises(TypeError):
            cf.srmr_to_rmsea("0.08", 4, 2)


# ---- test-equivalence.R -------------------------------------------------------------

def ref_rules(cfi, tli, rmsea, srmr, df, dd):
    c_, t_, r_, s_ = cfi > .95, tli > .95, rmsea < .06, srmr < .08
    s_sc = srmr < .06 * math.sqrt(df / dd)
    return {"cosmin": c_ or t_ or r_ or s_, "no_rmsea": c_ or t_ or s_,
            "srmr_df": c_ or t_ or r_ or s_sc,
            "two_of_three": (c_ + r_ + s_) >= 2, "all_three": c_ and r_ and s_,
            "hu_bentler": (cfi >= .95 and srmr <= .09) or (rmsea <= .06 and srmr <= .09)}


class TestEquivalence(unittest.TestCase):
    def test_reference(self):
        rng = random.Random(1)
        n = 5000
        X = [{"cfi": rng.uniform(.85, 1), "tli": rng.uniform(.85, 1),
              "rmsea": rng.uniform(0, .15), "srmr": rng.uniform(0, .15),
              "p": rng.randint(4, 20)} for _ in range(n)]
        for r in X:
            r["df"] = cf.cfa_df(r["p"])
        refs = [ref_rules(r["cfi"], r["tli"], r["rmsea"], r["srmr"], r["df"],
                          cf.n_moments(r["p"])) for r in X]
        for rule in ref_rules(.9, .9, .1, .1, 1, 1):
            got = cf.rate_cfa([r["cfi"] for r in X], [r["tli"] for r in X],
                              [r["rmsea"] for r in X], [r["srmr"] for r in X],
                              rule=rule, p=[r["p"] for r in X], df=[r["df"] for r in X])
            want = ["sufficient" if f[rule] else "insufficient" for f in refs]
            self.assertEqual(got, want, rule)
        self.assertEqual(len({f["cosmin"] for f in refs}), 2)
        self.assertEqual(len({f["all_three"] for f in refs}), 2)
        # with complete data both modes agree, and nothing is indeterminate
        a = cf.rerate(X, rules=RULES, missing="reported")
        b = cf.rerate(X, rules=RULES, missing="unknown")
        for rule in RULES:
            col = "rating_" + rule
            self.assertEqual([r[col] for r in a], [r[col] for r in b], rule)
            self.assertNotIn("indeterminate", [r[col] for r in a])


# ---- test-rate-cfa.R ----------------------------------------------------------------

def any_rep(*v):
    v = [x for x in v if x is not None]
    return None if not v else any(v)


def and_rep(*v):
    return None if any(x is None for x in v) else all(v)


REF_REPORTED = {
    "cosmin": lambda c, t, r, s: any_rep(c, t, r, s),
    "no_rmsea": lambda c, t, r, s: any_rep(c, t, s),
    "srmr_df": lambda c, t, r, s: any_rep(c, t, r, s),
    "two_of_three": lambda c, t, r, s: (
        None if sum(x is not None for x in (c, r, s)) < 2
        else sum(bool(x) for x in (c, r, s)) >= 2),
    "all_three": lambda c, t, r, s: and_rep(c, r, s),
    "hu_bentler": lambda c, t, r, s: any_rep(and_rep(c, s), and_rep(r, s)),
    "prinsen2016_a": lambda c, t, r, s: and_rep(any_rep(c, t), any_rep(r, s)),
    "prinsen2016_b": lambda c, t, r, s: any_rep(and_rep(any_rep(c, t), r), s),
}


def ref_unknown(f, c, t, r, s):
    x = [c, t, r, s]
    miss = [i for i, v in enumerate(x) if v is None]
    out = set()
    for combo in itertools.product((False, True), repeat=len(miss)):
        y = list(x)
        for i, v in zip(miss, combo):
            y[i] = v
        out.add(f(*y))
    return out.pop() if len(out) == 1 else None


def lab(v):
    return "indeterminate" if v is None else "sufficient" if v else "insufficient"


class TestRateCfa(WarnMixin, unittest.TestCase):
    def test_rule_names(self):
        self.assertEqual([r["rule"] for r in cf.cosmin_rules()], RULES)
        self.assertEqual(list(cf.cosmin_rules()[0]),
                         ["rule", "definition", "indices", "needs_p_df"])
        self.assertEqual([r["needs_p_df"] for r in cf.cosmin_rules()],
                         [False, False, True, False, False, False, False, False])

    def test_truth_table_81(self):
        vals = {"cfi": (0.97, 0.90), "tli": (0.97, 0.90),
                "rmsea": (0.04, 0.08), "srmr": (0.03, 0.10)}
        g = list(itertools.product((True, False, None), repeat=4))
        # R's expand.grid varies the first factor fastest
        g = [tuple(reversed(k)) for k in g]
        g.sort(key=lambda k: [[True, False, None].index(v) for v in reversed(k)])
        num = {name: [None if k[i] is None else vals[name][0 if k[i] else 1] for k in g]
               for i, name in enumerate(("cfi", "tli", "rmsea", "srmr"))}
        self.assertEqual(len(g), 81)
        for rule in RULES:
            for mode in ("reported", "unknown"):
                got = cf.rate_cfa(num["cfi"], num["tli"], num["rmsea"], num["srmr"],
                                  rule=rule, p=10, df=35, missing=mode)
                f = REF_REPORTED[rule]
                want = [lab(f(*k) if mode == "reported" else ref_unknown(f, *k)) for k in g]
                self.assertEqual(got, want, f"{rule} {mode}")

    def test_all_missing(self):
        for rule in RULES:
            for mode in ("reported", "unknown"):
                self.assertEqual(cf.rate_cfa(rule=rule, missing=mode), "indeterminate")

    def test_single_index(self):
        r = cf.rate_cfa
        self.assertEqual(r(cfi=0.97), "sufficient")
        self.assertEqual(r(cfi=0.90), "insufficient")
        self.assertEqual(r(cfi=0.90, missing="unknown"), "indeterminate")
        self.assertEqual(r(rmsea=0.05), "sufficient")
        self.assertEqual(r(rmsea=0.05, rule="no_rmsea"), "indeterminate")
        self.assertEqual(r(cfi=0.97, rule="two_of_three"), "indeterminate")
        self.assertEqual(r(tli=0.97, srmr=0.03, rule="two_of_three"), "indeterminate")
        self.assertEqual(r(cfi=0.97, rule="all_three"), "indeterminate")
        self.assertEqual(r(cfi=0.90, rule="all_three", missing="unknown"), "insufficient")
        self.assertEqual(r(cfi=0.97, rmsea=0.04, rule="two_of_three", missing="unknown"),
                         "sufficient")

    def test_equality_edges(self):
        r = cf.rate_cfa
        self.assertEqual(r(cfi=0.95), "insufficient")
        self.assertEqual(r(tli=0.95), "insufficient")
        self.assertEqual(r(rmsea=0.06), "insufficient")
        self.assertEqual(r(srmr=0.08), "insufficient")
        self.assertEqual(r(srmr=0.0799), "sufficient")
        self.assertEqual(r(cfi=0.95, srmr=0.09, rule="hu_bentler"), "sufficient")
        self.assertEqual(r(rmsea=0.06, srmr=0.09, rule="hu_bentler"), "sufficient")
        self.assertEqual(r(cfi=0.97, srmr=0.0901, rule="hu_bentler"), "insufficient")
        self.assertEqual(r(cfi=0.97, rule="hu_bentler"), "indeterminate")
        self.assertEqual(r(cfi=0.90, rmsea=0.04, srmr=0.09, rule="hu_bentler"), "sufficient")

    def test_srmr_df(self):
        r = cf.rate_cfa
        self.assertEqual(r(srmr=0.03, rule="srmr_df"), "indeterminate")
        self.assertEqual(r(srmr=0.03, rule="srmr_df", p=10), "indeterminate")
        self.assertEqual(r(srmr=0.03, rule="srmr_df", p=10, df=35), "sufficient")
        self.assertEqual(r(srmr=0.06, rule="srmr_df", p=10, df=35), "insufficient")
        self.assertEqual(r(srmr=0.06), "sufficient")
        self.assertEqual(r(srmr=0.03, rule="srmr_df", p=3, df=0), "indeterminate")
        self.assertEqual(r(cfi=0.97, rule="srmr_df"), "sufficient")

    def test_prinsen2016(self):
        r = cf.rate_cfa
        self.assertEqual(r(cfi=0.97, srmr=0.03, rule="prinsen2016_a"), "sufficient")
        self.assertEqual(r(cfi=0.97, rule="prinsen2016_a"), "indeterminate")
        self.assertEqual(r(cfi=0.90, srmr=0.03, rule="prinsen2016_a"), "insufficient")
        self.assertEqual(r(srmr=0.03, rule="prinsen2016_b"), "sufficient")
        self.assertEqual(r(cfi=0.97, srmr=0.10, rule="prinsen2016_b"), "insufficient")
        self.assertEqual(r(cfi=0.97, srmr=0.10, rule="prinsen2016_b", missing="unknown"),
                         "indeterminate")
        self.assertEqual(r(cfi=0.97, rule="prinsen2016_b"), "indeterminate")
        self.assertEqual(r(cfi=0.97, rmsea=0.04, srmr=0.10, rule="prinsen2016_b"),
                         "sufficient")
        self.assertEqual(r(cfi=0.97, rmsea=0.08, rule="prinsen2016_b"), "insufficient")
        self.assertEqual(r(cfi=0.97, rmsea=0.08, rule="prinsen2016_b", missing="unknown"),
                         "indeterminate")

    def test_input_handling(self):
        r = cf.rate_cfa
        self.assertEqual(r(cfi=[0.97, 0.90], srmr=0.10), ["sufficient", "insufficient"])
        self.assertEqual(r(cfi=[0.97, 0.90, 0.97, 0.90], srmr=[0.05, 0.10]),
                         ["sufficient", "insufficient", "sufficient", "insufficient"])
        self.assertEqual(r(cfi=[]), [])
        with self.assertRaises(ValueError):
            r(cfi=[0.9, 0.9], srmr=[0.1, 0.1, 0.1])
        self.assertEqual(r(cfi="n.r.", srmr=0.05), "sufficient")
        self.assertEqual(r(cfi="n.r."), "indeterminate")
        self.assertEqual(r(cfi="0.97"), "sufficient")
        self.assertEqual(r(cfi=" 0.97 "), "sufficient")
        self.assertEqual(r(cfi=""), "indeterminate")
        self.assertEqual(r(cfi=float("inf")), "indeterminate")
        self.assertEqual(r(cfi=float("nan")), "indeterminate")
        self.assertEqual(r(cfi=None), "indeterminate")
        self.assertEqual(r(cfi={"a": 0.97}), ["sufficient"])    # R named vector
        self.assertEqual(r(cfi=(0.97,)), ["sufficient"])
        self.assertEqual(r(cfi=True), "sufficient")             # as.numeric(TRUE) == 1
        self.assertEqual(r(cfi="1_0"), "indeterminate")         # R: NA
        self.assertEqual(r(cfi="0x1"), "sufficient")            # R: 1
        self.assertEqual(r(cfi="1e"), "sufficient")             # R: 1
        with self.assertRaises(ValueError):
            r(cfi=0.97, rule="nope")
        with self.assertRaises(ValueError):
            r(cfi=0.97, missing="nope")
        self.assertWarns_(r, cfi=1.3)
        self.assertWarns_(r, tli=1.25)
        self.assertWarns_(r, srmr=-0.01)
        self.assertWarns_(r, rmsea=-0.01)
        self.assertNoWarning(r, cfi=1.1, tli=1.15, rmsea=0, srmr=0)

    def test_numpy_pandas_input(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest("numpy not installed")
        self.assertEqual(cf.rate_cfa(cfi=np.array([0.97, np.nan])),
                         ["sufficient", "indeterminate"])
        self.assertEqual(cf.rate_cfa(cfi=np.float64(0.97)), "sufficient")
        try:
            import pandas as pd
        except ImportError:
            return
        s = pd.Series([0.97, None, 0.90], dtype="Float64")      # contains pd.NA
        self.assertEqual(cf.rate_cfa(cfi=s),
                         ["sufficient", "indeterminate", "insufficient"])


# ---- test-rerate.R ------------------------------------------------------------------

D = [dict(cfi=0.97, tli=None, rmsea=0.05, srmr=0.04, p=10, df=35),
     dict(cfi=0.91, tli=0.89, rmsea=0.05, srmr=0.09, p=6, df=9),
     dict(cfi=None, tli=None, rmsea=None, srmr=0.07, p=8, df=20),
     dict(cfi=None, tli=None, rmsea=0.08, srmr=None, p=5, df=5),
     dict(cfi=0.90, tli=None, rmsea=0.07, srmr=0.05, p=6, df=9)]


class TestRerate(unittest.TestCase):
    def check_out(self, out, col):
        return [r[col] for r in out]

    def test_list_of_dicts(self):
        out = cf.rerate(D)
        self.assertTrue(all(set(d) <= set(o) for d, o in zip(D, out)))
        for rule in ("cosmin", "no_rmsea", "srmr_df", "two_of_three", "all_three"):
            self.assertIn("rating_" + rule, out[0])
        self.assertEqual(self.check_out(out, "n_indices"), [3, 4, 1, 1, 3])
        self.assertEqual(self.check_out(out, "rmsea_rescue"), [False, True, None, False, False])
        self.assertEqual(self.check_out(out, "rating_cosmin"),
                         ["sufficient", "sufficient", "sufficient", "insufficient", "sufficient"])
        self.assertEqual(out[2]["rating_srmr_df"], "insufficient")
        self.assertEqual(self.check_out(out, "rating_all_three"),
                         ["sufficient", "insufficient", "indeterminate", "indeterminate",
                          "insufficient"])
        for rule in ("cosmin", "all_three"):
            self.assertEqual(self.check_out(out, "rating_" + rule),
                             cf.rate_cfa(*[[d[k] for d in D] for k in
                                           ("cfi", "tli", "rmsea", "srmr")],
                                         rule=rule, p=[d["p"] for d in D],
                                         df=[d["df"] for d in D]))
        self.assertNotIn("rating_cosmin", D[0])          # input not modified

    def test_custom_cols(self):
        d2 = [{"CFI": 0.97, "SRMR": 0.10}, {"CFI": 0.90, "SRMR": 0.10}]
        out2 = cf.rerate(d2, rules="cosmin", cols={"cfi": "CFI", "srmr": "SRMR"})
        self.assertEqual(self.check_out(out2, "rating_cosmin"), ["sufficient", "insufficient"])
        self.assertEqual(self.check_out(out2, "n_indices"), [2, 2])
        self.assertEqual(self.check_out(out2, "rmsea_rescue"), [None, None])
        out3 = cf.rerate(d2, rules=["cosmin"], missing="unknown",
                         cols={"cfi": "CFI", "srmr": "SRMR"})
        self.assertEqual(self.check_out(out3, "rating_cosmin"),
                         ["sufficient", "indeterminate"])
        with self.assertRaises(ValueError):
            cf.rerate(D, rules=["nope"])
        with self.assertRaises(ValueError):
            cf.rerate(D, missing="nope")
        with self.assertRaises(TypeError):
            cf.rerate({k: [d[k] for d in D] for k in D[0]})   # not a data frame
        self.assertEqual(cf.rerate([]), [])

    def test_pandas(self):
        try:
            import pandas as pd
        except ImportError:
            self.skipTest("pandas not installed")
        df = pd.DataFrame(D)
        out = cf.rerate(df)
        self.assertIsInstance(out, pd.DataFrame)
        self.assertNotIn("rating_cosmin", df.columns)    # input not modified
        self.assertEqual(out["n_indices"].tolist(), [3, 4, 1, 1, 3])
        self.assertEqual([None if v is pd.NA else v for v in out["rmsea_rescue"]],
                         [False, True, None, False, False])
        self.assertEqual(out["rating_all_three"].tolist(),
                         self.check_out(cf.rerate(D), "rating_all_three"))
        out2 = cf.rerate(pd.DataFrame({"CFI": [0.97, 0.90], "SRMR": [0.10, 0.10]}),
                         rules="cosmin", cols={"cfi": "CFI", "srmr": "SRMR"})
        self.assertEqual(out2["rating_cosmin"].tolist(), ["sufficient", "insufficient"])
        self.assertTrue(out2["rmsea_rescue"].isna().all())


class TestReadme(unittest.TestCase):
    def test_readme_example(self):
        self.assertEqual(len(cf.cosmin_rules()), 8)
        self.assertEqual(cf.rate_cfa(cfi=0.93, rmsea=0.05), "sufficient")
        # The R README comments this call as "indeterminate", but R's
        # rate_cfa() returns "insufficient" (two of three indices reported,
        # only one passes). "indeterminate" is the missing = "unknown" result.
        self.assertEqual(cf.rate_cfa(cfi=0.93, rmsea=0.05, rule="two_of_three"),
                         "insufficient")
        self.assertEqual(cf.rate_cfa(cfi=0.93, rmsea=0.05, rule="two_of_three",
                                     missing="unknown"), "indeterminate")
        self.assertEqual(cf.rate_cfa(cfi=0.93, missing="unknown"), "indeterminate")
        self.assertEqual(cf.__version__, "0.1.0")


if __name__ == "__main__":
    unittest.main()
