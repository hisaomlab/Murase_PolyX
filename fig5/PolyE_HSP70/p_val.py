import pandas as pd
import itertools
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# 読み込み
df = pd.read_excel("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig5/PolyE_HSP70/20250921_PolyE_HSP70/p_val.xlsx")

groups = ["Vector", "EGFP", "PolyE"]

pairs = []
pvals = []

for g1, g2 in itertools.combinations(groups, 2):
    x = df[g1].dropna()
    y = df[g2].dropna()
    pvals.append(ttest_ind(x, y, equal_var=False).pvalue)
    pairs.append(f"{g1} vs {g2}")

# Bonferroni
_, p_corr, _, _ = multipletests(pvals, method="fdr_bh")

for pair, p, pc in zip(pairs, pvals, p_corr):
    print(f"{pair}: p = {p:.5g}, p_adj = {pc:.5g}")