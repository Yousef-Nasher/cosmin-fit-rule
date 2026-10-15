"""Independent check of the Study 2 headline numbers.
Written in pandas and sharing no code with the R analysis; reads the frozen zip directly.
Run from study2/:  python3 code/verify_independent.py > results/verify_independent.txt"""
import io, zipfile
import numpy as np
import pandas as pd

with zipfile.ZipFile("data/study2_final_chunks.zip") as z:
    jobs = sorted(n for n in z.namelist() if n.startswith("final/") and n.endswith(".csv"))
    D = pd.concat([pd.read_csv(io.BytesIO(z.read(n))) for n in jobs], ignore_index=True)
    P = pd.read_csv(io.BytesIO(z.read("population.csv")))
P = P[P.est == "ULS"].set_index("cond")
D["cond"] = D.cell.str.replace(r"_n\d{4}$", "", regex=True)
D["conseq"] = D.cond.map(P.consequential).astype(bool)
assert len(jobs) == 1520 and len(D) == 152 * 500 * 3
assert D.groupby("cell").rep.nunique().eq(500).all()
D = D[(D.status == "ok") & (D.admiss == True)].copy()   # clean fits

w = D.arm == "wlsmv_ord"            # WLSMV: scaled indices; ML arms: ML indices
cfi = np.where(w, D["cfi.sc"], D.cfi)
tli = np.where(w, D["tli.sc"], D.tli)
rmsea = np.where(w, D["rmsea.sc"], D.rmsea)
c = np.nan_to_num(cfi) > .95
t = np.nan_to_num(tli) > .95
r = np.nan_to_num(rmsea, nan=9) < .06
s = D.srmr.fillna(9) < .08
D["cosmin"] = c | t | r | s
D["no_rmsea"] = c | t | s
for arm in ["ml_cont", "ml_ord_as_cont", "wlsmv_ord"]:
    for n in [100, 200, 500, 1000]:
        g = (D.arm == arm) & (D.n == n)
        false_fail = 100 * (~D.cosmin[g & (D.kind == "none")]).mean()
        false_pass = 100 * D.cosmin[g & D.conseq].mean()
        changed = int((D.cosmin != D.no_rmsea)[g].sum())
        print(f"{arm:15s} N={n:4d}  false-fail {false_fail:4.1f}  false-pass {false_pass:4.1f}  RMSEA-changes {changed}")
