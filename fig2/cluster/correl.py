import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["font.family"] = "Arial"

# ==== 1) データ読み込み ====
data = pd.read_excel(
    "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/data.xlsx",
    sheet_name="raw",
    index_col=0
)

# ==== 2) スピアマン相関 ====
corr = data.corr(method="spearman")
np.fill_diagonal(corr.values, 1.0)  # 対角を1に固定

# ==== 3) 上三角をマスク ====
mask = np.triu(np.ones_like(corr, dtype=bool), k=1)

plt.figure(dpi=300, figsize=(14, 10))
sns.set(font_scale=1.6)

ax = sns.heatmap(
    corr,
    mask=mask,
    vmin=-1, vmax=1,
    cmap="RdBu_r",
    square=True,
    cbar_kws={"shrink": 0.8},
    annot=True,      # 数値を表示
    fmt=".2f",       # 小数点以下2桁
    linewidths=0,    # セルの境界線なし
    linecolor=None,   # グリッド線なし
    annot_kws={"size": 10}
)

# 背景を白に
ax.set_facecolor("white")

# 軸まわりの体裁
ax.xaxis.tick_bottom()
ax.yaxis.tick_left()
plt.xticks(rotation=90, fontsize=15)
plt.yticks(rotation=0, fontsize=15)

plt.tight_layout()
plt.savefig(
    "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/cluster/NI_correlation.png",
    bbox_inches="tight",
    dpi=300
)