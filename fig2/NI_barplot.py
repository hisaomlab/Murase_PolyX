import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_excel("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/data.xlsx",sheet_name="raw", index_col=0)
#order = ['K', 'R', 'W', 'M', 'F', 'Y', 'V', 'I', 'L', 'C', 'P', 'T', 'H', 'G', 'A', 'D', 'E', 'S', 'N', 'Q']

# 平均と標準偏差を計算
mean_df = df.mean(axis=1)
std_df = df.std(axis=1)

log_df = np.log10(df)
var_df  = log_df.var(axis=1)

# CV[%] を計算（ゼロ除算回避）
cv_df = (std_df / mean_df) * 100
cv_df = cv_df.replace([np.inf, -np.inf], np.nan)

# DataFrameにまとめる
summary_df = pd.DataFrame({"mean": mean_df, "std": std_df, "var": var_df, "cv": cv_df})

# 平均値でソート（サンプル名も保持）
summary_df = summary_df.sort_values(by="mean", ascending=False)
sorted_labels = summary_df.index.tolist()  # ソート後の順序
#summary_df = summary_df.reindex(order)

# プロット準備
fig, ax = plt.subplots(figsize=(8, 3))

# NI 平均値と標準偏差の棒グラフ
bars = ax.bar(summary_df.index, summary_df["mean"], color=np.where(summary_df.index.isin(["Δ"]), "#7F7F7F", "#D0B9A4"), width=0.6, linewidth=0.5, edgecolor='black', alpha=0.7)
#ax.errorbar(summary_df.index, summary_df["mean"], yerr=summary_df["std"], fmt='none', ecolor='black', elinewidth=0.5, capsize=0)

# 生データのドットプロット
for i, label in enumerate(sorted_labels):
    ax.scatter([i]*df.shape[1], df.loc[label], color='black', s=10, edgecolor='black', linewidths=0.1, zorder=3, alpha=0.5)

# for i, label in enumerate(order):
#     if label in df.index: 
#         ax.scatter([i]*df.shape[1], df.loc[label], color='black', s=10, edgecolor='black', linewidths=0.1, zorder=3, alpha=0.5)

# Y軸の範囲設定
#ax.set_ylim(0, 300)
ax.set_ylim(1, 100000)
ax.set_yscale('log')

# 保存
plt.savefig("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/NI_plot.png", dpi=350)

summary_df.index.name = "sample"

# Excel
summary_path_xlsx = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/NI_summary.xlsx"
summary_df.to_excel(summary_path_xlsx, sheet_name="summary")