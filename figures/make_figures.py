"""Figures for Studies 2, 2b and 3 (PNG, 300 dpi).
Run from the repository root:  python3 figures/make_figures.py
Reads study2/results/table_rules_by_N.csv, study2b/data/study2b_chunks.zip and the
Study 3 results in study3/results/. Needs pandas, numpy and matplotlib."""
import io
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#888888",
                     "axes.labelcolor": "#333333", "xtick.color": "#444444", "ytick.color": "#444444"})
OUT = "figures/"
INK, GRID = "#222222", "#e6e6e3"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]


def save(fig, name):
    fig.savefig(OUT + name + ".png", dpi=300, bbox_inches="tight", pad_inches=0.05,
                metadata={"Software": None})
    plt.close(fig)


# ---- Study 2: error rates of the rule and three alternatives, WLSMV scaled indices ----
T = pd.read_csv("study2/results/table_rules_by_N.csv")
T = T[(T.version == "scaled") & (T.arm == "wlsmv_ord")]
rules = [("COSMIN", "COSMIN rule"), ("SRMR_df_scaled", "SRMR cutoff scaled by df"),
         ("two_of_three", "Two of three"), ("all_three", "All three")]
Ns = [100, 200, 500, 1000]
x = np.arange(4)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))
for ax, col, title in [(axes[0], "false_fail", "A. Correct models failed"),
                       (axes[1], "false_pass", "B. Consequential misspecification passed")]:
    for (rule, label), colour, marker in zip(rules, SERIES, ["o", "s", "^", "D"]):
        y = [T[(T.rule == rule) & (T.n == n)][col].iloc[0] for n in Ns]
        ax.plot(x, y, color=colour, lw=2, marker=marker, ms=5, label=label,
                markeredgecolor="white", markeredgewidth=1)
    ax.set_xticks(x); ax.set_xticklabels([str(n) for n in Ns]); ax.set_xlabel("Sample size (N)")
    ax.set_ylabel("Percent of fits"); ax.set_ylim(0, 80 if col == "false_pass" else 25)
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    ax.set_title(title, loc="left", fontsize=9, color=INK)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=8)
fig.tight_layout(rect=(0, 0.08, 1, 1))
save(fig, "study2_error_rates")

# ---- Study 2b: RMSEA rescue rate for correct one-factor models, WLSMV ------------------
with zipfile.ZipFile("study2b/data/study2b_chunks.zip") as z:
    jobs = sorted(n for n in z.namelist() if n.startswith("chunks/") and n.endswith(".csv"))
    D = pd.concat([pd.read_csv(io.BytesIO(z.read(n))) for n in jobs], ignore_index=True)
D = D[(D.status == "ok") & (D.admiss == True) & (D.arm == "wlsmv_ord") & (D.kind == "none")].copy()
c = D["cfi.sc"].fillna(0) > .95
t = D["tli.sc"].fillna(0) > .95
r = D["rmsea.sc"].fillna(9) < .06
s = D.srmr.fillna(9) < .08
D["rescue"] = 100 * (r & ~c & ~t & ~s)
formats = ["c2_sym", "c2_skew", "c3_sym", "c3_skew", "c5_sym", "c5_skew"]
labels = {"c2_sym": "2 categories, symmetric", "c2_skew": "2 categories, skewed",
          "c3_sym": "3 categories, symmetric", "c3_skew": "3 categories, skewed",
          "c5_sym": "5 categories, symmetric", "c5_skew": "5 categories, skewed"}
cmap = LinearSegmentedColormap.from_list("b", ["#f7f9fc", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0), sharey=True)
for i, load in enumerate([.45, .70]):
    for j, p in enumerate([5, 10, 20]):
        ax = axes[i, j]
        g = D[(D.load == load) & (D.p == p)]
        M = g.pivot_table(index="ct", columns="n", values="rescue", aggfunc="mean").reindex(formats)
        ax.imshow(M.values, cmap=cmap, vmin=0, vmax=60, aspect="auto")
        for a in range(M.shape[0]):
            for b in range(M.shape[1]):
                v = M.values[a, b]
                ax.text(b, a, f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                        color="white" if v > 30 else INK)
        ax.set_xticks(range(3)); ax.set_xticklabels(["100", "200", "500"]); ax.set_yticks(range(6))
        ax.set_yticklabels([labels[f] for f in formats], fontsize=7.5)
        ax.set_title(f"{p} items, loadings {load:.2f}", fontsize=8.5, color=INK)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks(np.arange(-.5, 3), minor=True); ax.set_yticks(np.arange(-.5, 6), minor=True)
        ax.grid(which="minor", color="white", linewidth=2); ax.tick_params(which="minor", length=0)
        if i == 1:
            ax.set_xlabel("Sample size (N)")
fig.tight_layout()
save(fig, "study2b_rmsea_rescue")

# ---- Study 3: flow of records ---------------------------------------------------------
flow = pd.read_csv("study3/results/flow_counts.csv")
n = dict(zip(flow.stage + "|" + flow.detail, flow.n))
ex = flow[(flow.detail == "full text") & flow.stage.str.startswith("excluded:")]
fig, ax = plt.subplots(figsize=(7.2, 7.0))
ax.set_xlim(0, 100); ax.set_ylim(11, 100); ax.axis("off")


def box(x, y, w, h, text, fs=7.6, ha="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1",
                                fc="white", ec="#666666", lw=0.8))
    tx = x + w / 2 if ha == "center" else x + 1.2
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs, color=INK, linespacing=1.35)


def arrow(x1, y1, x2, y2):
    ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", color="#666666", lw=0.8))


L, W, R, RW = 1, 45, 52, 47
box(L, 90, W, 9, f"Records identified\nPubMed: {n['records identified|PubMed search']:,}; "
                 f"COSMIN database: {n['records identified|COSMIN database export']:,}")
box(R, 90, RW, 9, f"Duplicates removed: {n['duplicates removed|both sources']}")
arrow(L + W, 94.5, R, 94.5)
box(L, 77, W, 9, f"Records screened on title and abstract\n(n = {n['records screened on title and abstract|both sources']:,})")
arrow(L + W / 2, 90, L + W / 2, 86.3)
ta = flow[flow.detail == "title and abstract"]
ta = ta[ta.stage.str.startswith("excluded")]
box(R, 75, RW, 13, f"Excluded: {ta.n.sum()}\n" + "\n".join(
    f"{s.replace('excluded: ', '')[0].upper() + s.replace('excluded: ', '')[1:]}: {k}" for s, k in zip(ta.stage, ta.n)),
    fs=7, ha="left")
arrow(L + W, 81.5, R, 81.5)
box(L, 63, W, 9, f"Sampling frame: {n['records in the sampling frame|frame']:,} reviews\n"
                 f"(random order within publication year)")
arrow(L + W / 2, 77, L + W / 2, 72.3)
box(R, 63, RW, 9, f"Removed, publication year before 2018\nfound from DOI or URL: "
                  f"{n['removed: publication year before 2018 found from DOI or URL|frame']}")
arrow(L + W, 67.5, R, 67.5)
box(L, 49, W, 9, f"Reached in draw order: {n['records reached in draw order|frame']}\n"
                 f"(100 drawn, then replacements)")
arrow(L + W / 2, 63, L + W / 2, 58.3)
box(R, 49, RW, 9, f"Not reached (year quota filled): {n['records not reached (stratum quota filled)|frame']}\n"
                   f"Full text not obtained: {n['full text not obtained|reached']}", fs=7, ha="left")
arrow(L + W, 53.5, R, 53.5)
box(L, 35, W, 9, f"Full texts assessed for eligibility\n(n = {n['full text assessed|reached']})")
arrow(L + W / 2, 49, L + W / 2, 44.3)
SHORT = {
    "no numeric CFA fit index in main text or supplementary files": "No fit value in main text or supplements",
    "no numeric CFA fit index in main text; supplementary files could not be obtained": "No fit value in main text,\n   supplements not obtained",
    "COSMIN criteria not applied to rate measurement properties": "COSMIN criteria not applied",
    "no confirmatory factor analysis or structural validity results": "No CFA or structural validity results",
    "structural validity assessed by exploratory analysis only": "Exploratory factor analysis only",
    "no numeric CFA fit index in main text; no supplementary files exist": "No fit value in main text,\n   no supplements exist",
    "published before 2018 (outside the frame)": "Published before 2018",
    "protocol, not a completed review": "Protocol"}
lines = [SHORT[s.replace("excluded: ", "")] for s in ex.stage]
box(R, 17, RW, 24, f"Excluded: {ex.n.sum()}\n" + "\n".join(
    f"{s}: {k}" for s, k in zip(lines, ex.n)), fs=7, ha="left")
arrow(L + W, 39.5, R, 39.5)
box(L, 20, W, 9, f"Reviews included: {n['reviews included|full text']}\n"
                 f"{n['study-level CFA results|included reviews']} study-level CFA results")
arrow(L + W / 2, 35, L + W / 2, 29.3)
save(fig, "study3_flow")

# ---- Study 3: indices reported and re-rating under the repaired rules -------------------
q2 = pd.read_csv("study3/results/q2_indices.csv")
ch = pd.read_csv("study3/results/q4_result_changes.csv")
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axes[0]
ax.bar(q2.n_indices.astype(str), q2.percent, color=SERIES[0], width=0.65)
for xi, (v, k) in enumerate(zip(q2.percent, q2.n)):
    ax.text(xi, v + 1, f"{k}", ha="center", fontsize=7.5, color=INK)
ax.set_xlabel("Indices reported (CFI, TLI, RMSEA, SRMR)"); ax.set_ylabel("Percent of CFA results")
ax.set_ylim(0, 50); ax.set_title("A. Indices reported per result", loc="left", fontsize=9, color=INK)
ax = axes[1]
names = {"no_rmsea": "RMSEA clause\ndropped", "srmr_df": "SRMR cut-off\nscaled by df", "two_of_three": "Two of\nthree"}
y = np.arange(len(ch))[::-1]
base = ch.n_sufficient_as_written
ins = 100 * ch.to_insufficient / base
ind = 100 * ch.to_indeterminate / base
ax.barh(y, ins, color=SERIES[1], height=0.55, label="Insufficient")
ax.barh(y, ind, left=ins, color="#bdbdbd", height=0.55, label="Indeterminate")
for yi, a, b in zip(y, ins, ind):
    ax.text(a + b + 1, yi, f"{a + b:.0f}%", va="center", fontsize=7.5, color=INK)
ax.set_yticks(y); ax.set_yticklabels([names[r] for r in ch.rule], fontsize=8)
ax.set_xlim(0, 75); ax.set_xlabel(f"Percent of the {base.iloc[0]} results sufficient under the rule as written")
ax.set_title("B. No longer sufficient under a repaired rule", loc="left", fontsize=9, color=INK)
ax.legend(frameon=False, fontsize=7.5, loc="center right")
for a in axes:
    a.grid(axis="y" if a is axes[0] else "x", color=GRID, lw=0.8); a.set_axisbelow(True)
    for side in ["top", "right"]:
        a.spines[side].set_visible(False)
fig.tight_layout()
save(fig, "study3_rerating")
print("done")
