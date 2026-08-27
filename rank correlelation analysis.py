"""
rank_correlation_analysis.py

Priority-2 investigation: do the FLAGGED diagnostics (UBF, Brown2, NCSU1,
-Ri, F2D) rank the same grid cells as extreme as the WORKING shear-based
ones (TI1, TI2, VWS, Endlich)?

If yes → the unit/magnitude flags from the W&J comparison are cosmetic;
Prosser's percentile-threshold framework doesn't care about absolute values.
If no → those diagnostics genuinely produce different physics, and we
need to handcode them before trusting them in the trend analysis.

Analysis:
  1. Spearman rank correlation matrix across all 21 diagnostics (day mean).
  2. Top-5% and top-1% Jaccard overlap — Prosser's actual test.
  3. Time stability: does the correlation hold across all 8 time steps?
  4. Hierarchical clustering: which diagnostics form natural groups?
  5. Verdict per flagged diagnostic.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform


OUT = Path("/home/claude/work/priority2_out")
OUT.mkdir(parents=True, exist_ok=True)

ds = xr.open_dataset("/home/claude/work/diagnostics.nc")
names = list(ds.data_vars)
N = len(names)

# ---------------------------------------------------------------------------
# 1. Build a (N × n_grid_points) matrix of |day-mean| values
# ---------------------------------------------------------------------------
flat = np.stack([np.abs(ds[k].mean("time")).values.ravel() for k in names])
# Drop non-finite columns (edge NaNs from gradients)
mask = np.isfinite(flat).all(axis=0)
flat = flat[:, mask]
print(f"Analysis matrix: {N} diagnostics × {flat.shape[1]} grid cells")

# ---------------------------------------------------------------------------
# 2. Spearman rank correlation matrix
# ---------------------------------------------------------------------------
print("\n>>> Computing Spearman rank correlations (day mean)...")
rho = np.zeros((N, N))
for i in range(N):
    for j in range(i, N):
        r, _ = spearmanr(flat[i], flat[j])
        rho[i, j] = rho[j, i] = r

rho_df = pd.DataFrame(rho, index=names, columns=names)
rho_df.to_csv(OUT / "spearman_daymean.csv")

# ---------------------------------------------------------------------------
# 3. Top-percentile Jaccard overlap
#    (Prosser's method: flag cells above threshold, count overlap.)
# ---------------------------------------------------------------------------
def jaccard(a, b, q):
    ta = np.quantile(a, q)
    tb = np.quantile(b, q)
    A = a >= ta
    B = b >= tb
    return (A & B).sum() / max((A | B).sum(), 1)

print(">>> Computing top-p95 and top-p99 Jaccard overlaps...")
J95 = np.zeros((N, N))
J99 = np.zeros((N, N))
for i in range(N):
    for j in range(i, N):
        J95[i, j] = J95[j, i] = jaccard(flat[i], flat[j], 0.95)
        J99[i, j] = J99[j, i] = jaccard(flat[i], flat[j], 0.99)

J95_df = pd.DataFrame(J95, index=names, columns=names)
J99_df = pd.DataFrame(J99, index=names, columns=names)
J95_df.to_csv(OUT / "jaccard_p95.csv")
J99_df.to_csv(OUT / "jaccard_p99.csv")

# ---------------------------------------------------------------------------
# 4. Time-stability: how consistent is the correlation across 8 timesteps?
# ---------------------------------------------------------------------------
print(">>> Time-stability across the 8 timesteps...")
n_time = ds.sizes["time"]
rho_time = np.zeros((N, N, n_time))
for t in range(n_time):
    flat_t = np.stack([np.abs(ds[k].isel(time=t)).values.ravel() for k in names])
    mask_t = np.isfinite(flat_t).all(axis=0)
    flat_t = flat_t[:, mask_t]
    for i in range(N):
        for j in range(i, N):
            r, _ = spearmanr(flat_t[i], flat_t[j])
            rho_time[i, j, t] = rho_time[j, i, t] = r

rho_time_std = rho_time.std(axis=2)  # dispersion across time
rho_time_min = rho_time.min(axis=2)  # worst timestep

# ---------------------------------------------------------------------------
# 5. Verdict tables — focus on flagged vs working
# ---------------------------------------------------------------------------
FLAGGED  = ["negative_richardson", "ubf", "brown2", "f2d", "ncsu1"]
WORKING_SHEAR = ["vertical_wind_shear", "endlich", "ti1", "ti2", "ngm1"]

def verdict_row(flag_name):
    i = names.index(flag_name)
    rows = []
    for w in WORKING_SHEAR:
        j = names.index(w)
        rows.append({
            "flagged": flag_name,
            "working_shear": w,
            "Spearman ρ (day mean)": rho[i, j],
            "Spearman ρ min across 8 timesteps": rho_time_min[i, j],
            "Jaccard @ top-5%": J95[i, j],
            "Jaccard @ top-1%": J99[i, j],
        })
    return rows

all_rows = []
for f in FLAGGED:
    all_rows.extend(verdict_row(f))
verdict_df = pd.DataFrame(all_rows)
verdict_df.to_csv(OUT / "flagged_vs_working_shear.csv", index=False)

print("\n>>> Verdict — flagged diagnostics vs working shear diagnostics")
pd.set_option("display.max_colwidth", None)
pd.set_option("display.width", 200)
pd.set_option("display.float_format", lambda x: f"{x:+.3f}")
print(verdict_df.to_string(index=False))

# ---------------------------------------------------------------------------
# 6. Also check flagged vs each other (do they at least agree with each
#    other, even if not with the shear cluster?)
# ---------------------------------------------------------------------------
print("\n>>> Flagged vs flagged (Spearman ρ, day mean)")
flag_rho = pd.DataFrame(index=FLAGGED, columns=FLAGGED, dtype=float)
for a in FLAGGED:
    for b in FLAGGED:
        flag_rho.loc[a, b] = rho[names.index(a), names.index(b)]
print(flag_rho.to_string())

# ---------------------------------------------------------------------------
# 7. Hierarchical clustering
# ---------------------------------------------------------------------------
# Convert to distance: 1 - |ρ|
D = 1 - np.abs(rho)
np.fill_diagonal(D, 0)
D = (D + D.T) / 2  # enforce exact symmetry
D_condensed = squareform(D, checks=False)
Z = linkage(D_condensed, method="average")

# Cluster into 5 groups
clusters = fcluster(Z, t=5, criterion="maxclust")
cluster_df = pd.DataFrame({"diagnostic": names, "cluster": clusters}).sort_values("cluster")
cluster_df.to_csv(OUT / "clusters.csv", index=False)
print("\n>>> Hierarchical clusters (5 groups, 1-|ρ| distance, average linkage)")
for c in sorted(cluster_df["cluster"].unique()):
    members = cluster_df[cluster_df["cluster"] == c]["diagnostic"].tolist()
    print(f"  Cluster {c}: {', '.join(members)}")

# ---------------------------------------------------------------------------
# 8. Plots
# ---------------------------------------------------------------------------
# Reorder for visual clustering
order = np.argsort(clusters)
names_ord = [names[i] for i in order]
rho_ord = rho[np.ix_(order, order)]
J95_ord = J95[np.ix_(order, order)]

fig, axes = plt.subplots(1, 2, figsize=(18, 8))
im1 = axes[0].imshow(rho_ord, cmap="RdBu_r", vmin=-1, vmax=1)
axes[0].set_title("Spearman rank correlation ρ (day mean, |diag|)\nordered by hierarchical clusters", fontsize=11)
axes[0].set_xticks(range(N)); axes[0].set_yticks(range(N))
axes[0].set_xticklabels(names_ord, rotation=90, fontsize=8)
axes[0].set_yticklabels(names_ord, fontsize=8)
plt.colorbar(im1, ax=axes[0], fraction=0.04)

im2 = axes[1].imshow(J95_ord, cmap="magma", vmin=0, vmax=1)
axes[1].set_title("Jaccard overlap @ top-5% cells\n(Prosser's exceedance framework)", fontsize=11)
axes[1].set_xticks(range(N)); axes[1].set_yticks(range(N))
axes[1].set_xticklabels(names_ord, rotation=90, fontsize=8)
axes[1].set_yticklabels(names_ord, fontsize=8)
plt.colorbar(im2, ax=axes[1], fraction=0.04)
plt.tight_layout()
plt.savefig(OUT / "correlation_matrices.png", dpi=140, bbox_inches="tight")
plt.close()
print(f"\nSaved: {OUT / 'correlation_matrices.png'}")

# Verdict-focused subplot: flagged (rows) vs working shear (cols)
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
sub_rho = np.array([[rho[names.index(f), names.index(w)] for w in WORKING_SHEAR] for f in FLAGGED])
sub_J   = np.array([[J95[names.index(f), names.index(w)] for w in WORKING_SHEAR] for f in FLAGGED])
im1 = axes[0].imshow(sub_rho, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axes[0].set_xticks(range(len(WORKING_SHEAR))); axes[0].set_xticklabels(WORKING_SHEAR, rotation=45, ha="right", fontsize=9)
axes[0].set_yticks(range(len(FLAGGED))); axes[0].set_yticklabels(FLAGGED, fontsize=9)
axes[0].set_title("Spearman ρ:  flagged (rows)  ×  working shear (cols)", fontsize=10)
for i in range(len(FLAGGED)):
    for j in range(len(WORKING_SHEAR)):
        axes[0].text(j, i, f"{sub_rho[i,j]:+.2f}", ha="center", va="center",
                     fontsize=8, color="white" if abs(sub_rho[i,j]) > 0.5 else "black")
plt.colorbar(im1, ax=axes[0], fraction=0.04)

im2 = axes[1].imshow(sub_J, cmap="magma", vmin=0, vmax=1, aspect="auto")
axes[1].set_xticks(range(len(WORKING_SHEAR))); axes[1].set_xticklabels(WORKING_SHEAR, rotation=45, ha="right", fontsize=9)
axes[1].set_yticks(range(len(FLAGGED))); axes[1].set_yticklabels(FLAGGED, fontsize=9)
axes[1].set_title("Jaccard @ top-5%:  flagged  ×  working shear", fontsize=10)
for i in range(len(FLAGGED)):
    for j in range(len(WORKING_SHEAR)):
        axes[1].text(j, i, f"{sub_J[i,j]:.2f}", ha="center", va="center",
                     fontsize=8, color="white" if sub_J[i,j] > 0.5 else "black")
plt.colorbar(im2, ax=axes[1], fraction=0.04)
plt.tight_layout()
plt.savefig(OUT / "flagged_vs_shear.png", dpi=140, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT / 'flagged_vs_shear.png'}")

print("\nAll outputs in:", OUT)