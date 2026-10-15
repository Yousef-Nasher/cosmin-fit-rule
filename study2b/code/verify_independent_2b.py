"""Independent check of the Study 2b expectation counts and key cells.
Written in pandas and sharing no code with the R analysis; reads the frozen zip directly.
Run from study2b/:  python3 code/verify_independent_2b.py > results/verify_independent_2b.txt"""
import io, zipfile
import pandas as pd

with zipfile.ZipFile("data/study2b_chunks.zip") as z:
    jobs = sorted(n for n in z.namelist() if n.startswith("chunks/") and n.endswith(".csv"))
    D = pd.concat([pd.read_csv(io.BytesIO(z.read(n))) for n in jobs], ignore_index=True)
assert len(jobs) == 1152 and len(D) == 288 * 200 * 2
assert D.groupby("cell").rep.nunique().eq(200).all()
# clean WLSMV fits of correct models, scaled indices
D = D[(D.status == "ok") & (D.admiss == True) & (D.arm == "wlsmv_ord") & (D.kind == "none")].copy()
c = D["cfi.sc"].fillna(0) > .95
t = D["tli.sc"].fillna(0) > .95
r = D["rmsea.sc"].fillna(9) < .06
s = D.srmr.fillna(9) < .08
D["srmr_fail"] = ~s
D["rescue"] = r & ~c & ~t & ~s          # RMSEA is the only clause that passes
g = D.groupby(["ct", "p", "load", "n"])[["srmr_fail", "rescue"]].mean().reset_index()
e1 = g[g.ct.isin(["c2_sym", "c2_skew", "c3_skew"]) & (g.n <= 200)]
print("E1", int((e1.srmr_fail > .2).sum()), "of", len(e1),
      " E2", int((e1.rescue > .05).sum()), "of", len(e1), " max rescue", round(100 * e1.rescue.max(), 1))
e3 = g[g.ct == "c5_sym"]
print("E3", int((e3.rescue < .01).sum()), "of", len(e3), " max", round(100 * e3.rescue.max(), 1))
for row in g[(g.ct == "c2_skew") & (g.p == 20) & (g.load == .45)].itertuples():
    print(f"c2_skew p=20 load=0.45 N={row.n:3d}  SRMR fails {100 * row.srmr_fail:5.1f}  RMSEA rescue {100 * row.rescue:5.1f}")
