import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
plt.rcParams["font.family"] = "Arial"
from ipywidgets import interact, FloatSlider

#毒性データ取得
all_df = pd.read_excel("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/NI/NI.xlsx",sheet_name="All_NI_Data", index_col=0)
all_df = all_df.drop("Δ", axis=0)
all_df = all_df.T

# 赤青のグラデーションを作成
custom_colors = [(0, "#012346"), (0.5, "white"), (1, "#FBD607")]
C = mcolors.LinearSegmentedColormap.from_list("custom_cmap", custom_colors)

#毒性クラスタリング
AA_loc = ["#f3e55b","#404040","#4cb88c","#4cb88c","#f3e55b","#404040","#ea5f65","#f3e55b","#ea5f65","#f3e55b",
                   "#f3e55b","#5387c5","#404040","#5387c5","#ea5f65","#5387c5","#5387c5", "#f3e55b","#f3e55b","#f3e55b"]

ax = sns.clustermap(all_df,method='average', metric='euclidean',
                                    linewidths=0,
                                    cmap=C,
                                    center=10000,
                                    vmin=0,        # ここで最小値を 0 に設定
                                    vmax=20000,         # ここで最大値を 2 に設定
                                    #row_colors=AA_loc,
                                    #col_colors=AA_loc,
                                    dendrogram_ratio=(0.15, 0.15),
                                    cbar_pos=(1.1, 0.3, 0.03, 0.4),
                                    figsize=(8, 3))

# 軸のラベルのテキストを指定する（文字の大きさを縦横軸で変えるため）
yticklabels = ax.ax_heatmap.get_yticklabels()
ax.ax_heatmap.set_yticklabels(yticklabels, fontsize=14)

xticklabels = ax.ax_heatmap.get_xticklabels()
ax.ax_heatmap.set_xticklabels(xticklabels, fontsize=12)

plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=0)
plt.setp(ax.ax_heatmap.get_yticklabels(), rotation=0)

# デンドログラムを非表示にする
ax.ax_row_dendrogram.set_visible(False)
ax.ax_col_dendrogram.set_visible(False)


plt.savefig("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/NI/cluster.png",bbox_inches='tight', dpi = 350)