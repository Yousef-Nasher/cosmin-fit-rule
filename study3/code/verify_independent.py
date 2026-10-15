"""Independent check of the Study 3 headline numbers.
Written in pandas and sharing no code with the R analysis or the cosminfit package.
Run from study3/:  python3 code/verify_independent.py > results/verify_independent.txt"""
import numpy as np
import pandas as pd

rev = pd.read_csv("data/study3_reviews.csv")
cfa = pd.read_csv("data/study3_cfa_results.csv", keep_default_na=False, na_values=[""])
ins = pd.read_csv("data/study3_instruments.csv", keep_default_na=False, na_values=[""])
for c in ["value_flags", "rating", "rating_scope", "instrument_group"]:
    cfa[c] = cfa[c].fillna("")
for c in ["sv", "ic", "recommendation", "sv_quality"]:
    ins[c] = ins[c].fillna("")

print("Q1 wording:", rev.wording.value_counts().reindex(
    ["disjunctive", "conjunctive", "mixed", "availability", "not stated"], fill_value=0).to_dict())

r = cfa[(cfa.level == "study") & ~cfa["value_flags"].str.contains("range")].copy()
IDX = ["cfi", "tli", "rmsea", "srmr"]
r["k"] = r[IDX].notna().sum(axis=1)
print("CFA results:", len(r), "reviews:", r.review_id.nunique())
print("Q2 indices per result:", r.k.value_counts().sort_index().to_dict(), "median", r.k.median())

# pass (1), fail (0) or not reported (NaN) for each clause
P = pd.DataFrame({
    "cfi": np.where(r.cfi.isna(), np.nan, r.cfi > .95),
    "tli": np.where(r.tli.isna(), np.nan, r.tli > .95),
    "rmsea": np.where(r.rmsea.isna(), np.nan, r.rmsea < .06),
    "srmr": np.where(r.srmr.isna(), np.nan, r.srmr < .08),
}, index=r.index)


def any_of(cols):
    """OR over the reported clauses: sufficient / insufficient / indeterminate."""
    x = P[cols]
    out = np.where(x.eq(1).any(axis=1), "sufficient", np.where(x.notna().any(axis=1), "insufficient", "indeterminate"))
    return pd.Series(out, index=r.index)


# df-scaled SRMR clause: reported df, otherwise one-factor df from the item count
p = pd.to_numeric(r["items"], errors="coerce")
df = r["df"].astype(float).where(r["df"].notna(), np.where(p >= 4, p * (p - 3) / 2, np.nan))
cut = .06 * np.sqrt(df / (p * (p + 1) / 2))
srmr_df = np.where(r.srmr.isna(), np.nan, np.where(cut.isna(), np.nan, r.srmr < cut))
srmr_unknown = r.srmr.notna() & (r.srmr < .08) & cut.isna()

R = pd.DataFrame(index=r.index)
R["cosmin"] = any_of(["cfi", "tli", "rmsea", "srmr"])
R["no_rmsea"] = any_of(["cfi", "tli", "srmr"])
Q = P[["cfi", "tli", "rmsea"]].copy()
Q["s"] = srmr_df
pass_any = Q.eq(1).any(axis=1)
eval_any = Q.notna().any(axis=1)
R["srmr_df"] = np.where(pass_any, "sufficient",
                        np.where(srmr_unknown, "indeterminate",
                                 np.where(eval_any, "insufficient", "indeterminate")))
t = P[["cfi", "rmsea", "srmr"]]
R["two_of_three"] = np.where(t.eq(1).sum(axis=1) >= 2, "sufficient",
                             np.where(t.notna().sum(axis=1) >= 2, "insufficient", "indeterminate"))
for rule in R:
    print(f"Q4 {rule:13s}", R[rule].value_counts().reindex(["sufficient", "insufficient", "indeterminate"],
                                                         fill_value=0).to_dict())
suff = R.cosmin == "sufficient"
for rule in ["no_rmsea", "srmr_df", "two_of_three"]:
    print(f"Q4 {rule:13s} sufficient -> insufficient {int((R[rule][suff] == 'insufficient').sum())}, "
          f"-> not sufficient {int((R[rule][suff] != 'sufficient').sum())} of {int(suff.sum())}")

own = r.rating.isin(["+", "-"]) & (r.rating_scope == "result")
print("Printed + with no passing index:", int((own & (r.rating == "+") & (R.cosmin == "insufficient")).sum()),
      "| printed - with a passing index:", int((own & (r.rating == "-") & (R.cosmin == "sufficient")).sum()))

fail = (r.cfi <= .95) | (r.tli <= .95) | (r.rmsea >= .06) | (r.srmr >= .08)
q3 = (r.rating == "+") & (r.rating_scope == "result") & (r.k >= 2)
print(f"Q3 {int(fail[q3].sum())}/{int(q3.sum())} = {100 * fail[q3].mean():.1f}% from {r.review_id[q3].nunique()} reviews")

# instrument level: + when >= 75% of evaluable results are sufficient, - when <= 25%, +/- otherwise
r["group"] = np.where(r.instrument_group == "", r.instrument, r.instrument_group)


def pool(s):
    e = s[s != "indeterminate"]
    if len(e) == 0:
        return "?"
    f = (e == "sufficient").mean()
    return "+" if f >= .75 else "-" if f <= .25 else "+/-"


G = R.assign(review_id=r.review_id, group=r.group).groupby(["review_id", "group"]).agg(pool)
G = G.join(rev.set_index("review_id")[["cosmin_version", "recommendation_scheme"]], on="review_id")
G = G.join(ins.set_index(["review_id", "instrument_group"])[["ic", "recommendation", "sv_quality"]].rename_axis(
    ["review_id", "group"]))
for c in ["ic", "recommendation", "sv_quality"]:
    G[c] = G[c].fillna("")
print("Instruments:", len(G), "| + as written:", int((G.cosmin == "+").sum()))
abc = (G.cosmin_version == 1) & (G.recommendation_scheme == "A/B/C")
n_rec = int((abc & (G.recommendation != "")).sum())
for rule in ["no_rmsea", "srmr_df", "two_of_three"]:
    lost = (G.cosmin == "+") & (G[rule] != "+")
    ic = G.ic.isin(["+", "-"]) & lost
    a_b = abc & (G.recommendation == "A") & ic
    b_c = abc & (G.recommendation == "B") & (G.cosmin != "-") & (G[rule] == "-") & (G.sv_quality == "high")
    print(f"{rule:13s} instruments + -> not + {int(lost.sum())}; IC -> ? {int(ic.sum())}; "
          f"A -> B {int(a_b.sum())}; B -> C {int(b_c.sum())}; recommendations changed {int((a_b | b_c).sum())}/{n_rec}")
