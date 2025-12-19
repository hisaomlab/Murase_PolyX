import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

plt.rcParams["font.family"] = "Arial"

# 毒性データ取得
all_df = pd.read_excel(
    "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/NI_correlation_results.xlsx",
    sheet_name="Spearman",
    index_col=0
)

# 赤青のグラデーションを作成
custom_colors = [(0, "#012346"), (0.5, "white"), (1, "#FB3407")]
C = mcolors.LinearSegmentedColormap.from_list("custom_cmap", custom_colors)

# 毒性クラスタリング
ax = sns.clustermap(
    all_df,
    method='average',
    metric='euclidean',
    linewidths=0,
    cmap=C,
    center=0,
    vmin=-1,
    vmax=1,
    figsize=(5, 8)
)

# 軸のフォントサイズ・回転
ax.ax_heatmap.tick_params(axis='x', labelsize=8, rotation=30)
ax.ax_heatmap.tick_params(axis='y', labelsize=8, rotation=0)

# 保存
ax.savefig("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/cluster.png", bbox_inches='tight', dpi=350)