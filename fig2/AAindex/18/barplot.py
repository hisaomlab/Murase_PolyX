import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# データ読み込み
df = pd.read_excel('/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/NI_correlation_results.xlsx', sheet_name='Spearman')
cols = [c for c in df.columns if c != "Property"]

for col in cols:
    sub_melted = df.melt(id_vars='Property', value_vars=col, var_name='Method', value_name='Correlation')
    order = sub_melted.sort_values("Correlation", ascending=False)["Property"]
    
    plt.figure(figsize=(8, max(4, 0.4 * len(sub_melted))))

    ax = sns.barplot(
                    data=sub_melted,
                    y='Property',
                    x='Correlation',
                    color="#000000",
                    order=order
    )

    ax.set_title(col, fontsize=10, loc='center')
    ax.set_xlabel("")
    ax.set_ylabel("")

    # y軸を右側に移動
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

    plt.subplots_adjust(left=0.75, right=0.95, top=0.85, bottom=0.05)

    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top') 

    for d in ['right', 'left', 'top', 'bottom']:
        ax.spines[d].set_color('black')

    ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    ax.xaxis.grid(True, linestyle='--', color='lightgray')
    ax.yaxis.grid(False)

    plt.xlim(-1.0, 1.0)
    plt.savefig(f"/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/barplot/barplot_{col}.png",bbox_inches='tight', dpi=350)

#=========

# mean, std を Property ごとに計算
mean_df = df[cols].mean(axis=1)
std_df  = df[cols].std(axis=1)

# Property を含めた DataFrame を作る
summary_df = pd.DataFrame({"Property": df["Property"], "mean": mean_df, "std": std_df})
order = summary_df.sort_values("mean", ascending=False)["Property"]

plt.figure(figsize=(8, max(4, 0.4 * len(summary_df))))

ax = sns.barplot(data=summary_df, y='Property', x='mean', color="#CCCCCC", order=order)

# --- 各メソッドの値を散布（x=値, y=カテゴリ位置）---
for i, prop in enumerate(order):
    vals = df.loc[df["Property"] == prop, cols].values.flatten()
    ax.scatter(vals, np.full_like(vals, i, dtype=float), s=12, edgecolor='#000000', linewidths=0.1, zorder=4, alpha=0.5, color='#000000')

ax.set_xlabel("")
ax.set_ylabel("")

# y軸を右側に移動
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

plt.subplots_adjust(left=0.75, right=0.95, top=0.85, bottom=0.05)

ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top') 

for d in ['right', 'left', 'top', 'bottom']:
    ax.spines[d].set_color('black')

ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax.xaxis.grid(True, linestyle='--', color='lightgray')
ax.yaxis.grid(False)

plt.xlim(-1.0, 1.0)

plt.savefig("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/barplot_summary.png", bbox_inches='tight', dpi=350)