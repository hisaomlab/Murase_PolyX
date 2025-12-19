#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

# Ensure Matplotlib has a writable cache directory before importing Matplotlib.
MPL_CACHE_DIR = Path(__file__).with_name(".matplotlib_cache")
MPL_CACHE_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))

import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 15
plt.rcParams["lines.linewidth"] = 2.5    # 線の太さ
plt.rcParams["axes.linewidth"] = 1.5     # 枠線の太さ
import pandas as pd

INPUT_PATH = Path(__file__).with_name("6FP_NI_pangenom.xlsx")
OUTPUT_PATH = Path(__file__).with_name("scatter_NI_vs_PanGenome.png")
XLIM = (0.0005, 100)  # e.g., (1, 3e4)
YLIM = (0.01, 10e4)  # e.g., (1, 5e4)

df = pd.read_excel(INPUT_PATH)

# before = len(df)
# df = df[(df["Pan_genom_num_poly10X_"] > 0) & (df["NI"] > 0)]
# dropped = before - len(df)
# if dropped:
#     print(f"Dropped {dropped} rows with non-positive values for log scale.")

amino_acids = sorted(df["AminoAcid"].dropna().unique())
# cmap = plt.get_cmap("tab20", len(amino_acids) or 1)
# color_map = {aa: cmap(i) for i, aa in enumerate(amino_acids)}

aa_colors = {
    'A': "#f3b45b", 'C': '#404040', 'D': '#4cb88c', 'E': '#4cb88c',
    'F': '#f3b45b', 'G': '#404040', 'H': '#ea5f65', 'I': '#f3b45b',
    'K': '#ea5f65', 'L': '#f3b45b', 'M': '#f3b45b', 'N': '#5387c5',
    'P': '#404040', 'Q': "#5387c5", 'R': '#ea5f65', 'S': '#5387c5',
    'T': '#5387c5', 'V': '#f3b45b', 'W': '#f3b45b', 'Y': '#f3b45b'
}

fig, ax = plt.subplots(figsize=(6, 6))

for aa, group in df.groupby("AminoAcid"):
    ax.scatter(
        group["Pan_Num_Poly10X_per_strain"],
        group["6FPs_Relative neutrality"],
        label=aa,
        #color=color_map.get(aa, "gray"),
        color=aa_colors.get(aa, "gray"),  # ← ここを変更
        edgecolors="none",
        linewidths=0.5,
        alpha = 0.8,
        s=80,
    )

ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_xlabel("Pan_genom_num_poly10X")
#ax.set_ylabel("NI")
#ax.set_title("Pan genome vs NI by Amino Acid")

# 枠線を整理
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

if XLIM:
    ax.set_xlim(XLIM)
if YLIM:
    ax.set_ylim(YLIM)

#ax.legend(title="AminoAcid", bbox_to_anchor=(1.05, 1), loc="upper left")
fig.tight_layout()

fig.savefig(OUTPUT_PATH, dpi=350)
print(f"Saved scatter plot to {OUTPUT_PATH}")
