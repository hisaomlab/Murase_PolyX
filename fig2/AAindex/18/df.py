import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import math
plt.rcParams["font.family"] = "Arial"

colors = {"EGFP_TDH3": "#00FF22",
          "EGFP_WTC": "#00FF22",
          "moxGFP": "#00FF22",
          "mNeonGreen": "#00FF22",
          "Gamillus":    "#00FF22",
          "mScarlet-I":    "#FF00AE",
          "mCherry":     "#FF00AE"}

# 毒性データ取得
all_df = pd.read_excel("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/AAindex1_property_18.xlsx", sheet_name="AAindex1_property_18")

target_cols = ["EGFP","moxGFP","mNeonGreen","Gamillus","mScarlet-I","mCherry"]
id_cols = [c for c in all_df.columns if c not in target_cols]

df = all_df.melt(id_vars=id_cols, value_vars=target_cols, var_name="FP", value_name="NI")

output_path = '/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/NI_property.xlsx'
df.to_excel(output_path, index=False)



y_cols = [c for c in id_cols if c != "AminoAcid"]
n = len(y_cols)
fig, axes = plt.subplots(7, 4, figsize=(8,14))
axes = axes.flatten()

for ax, y_col in zip(axes, y_cols):
    for fp in target_cols:
        subset = df[df["FP"] == fp]
        ax.scatter(subset["NI"], subset[y_col], label=fp, 
                   #color=colors[fp], 
                   color="#000000", 
                   s=8, alpha=0.3, edgecolor="none")
    ax.set_title(y_col, fontsize=4)
    ax.set_xlim(0.01, 100000)
    ax.set_xscale('log')



# 余った subplot を削除（y_cols が nrows*ncols より少ない場合）
for ax in axes[len(y_cols):]:
    fig.delaxes(ax)

# 凡例をまとめて右外に
handles, labels = axes[0].get_legend_handles_labels()


plt.tight_layout()
plt.savefig('/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/NI_property.png', dpi=350)