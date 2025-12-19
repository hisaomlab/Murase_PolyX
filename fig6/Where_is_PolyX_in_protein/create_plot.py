import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === 色分けルール ===
def get_color_map(amino_acids):
    color_map = {}
    for aa in amino_acids:
        if aa in "QNST":
            color_map[aa] = "#5285C1"  # 青
        elif aa in "ED":
            color_map[aa] = "#4CB88C"  # 緑
        elif aa in "RHK":
            color_map[aa] = "#EB5F65"  # 赤
        elif aa in "CPG":
            color_map[aa] = "#404040"  # 黒
        else:
            color_map[aa] = "#F4E55B"  # 黄
    return color_map

# === データ読み込み ===
csv_path = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig6/Where_is_PolyX_in_protein/Where_is_polyX_results.csv"
df = pd.read_csv(csv_path)
df = df.dropna(subset=["Avg pLDDT"])

amino_acids = sorted(df["PolyX Type"].unique())
palette = get_color_map(amino_acids)

# === Figure とサブプロット作成（4行×5列） ===
fig, axes = plt.subplots(4, 6, figsize=(20, 10))
axes = axes.flatten()  # 2D → 1D に変換

for i, aa in enumerate(amino_acids):
    ax = axes[i]
    subset = df[df["PolyX Type"] == aa]
    
    sns.scatterplot(
        data=subset,
        x="Length",
        y="Avg pLDDT",
        color="#000000",
        #edgecolor=palette.get(aa, "#000000"),
        alpha=0.5,
        s=20,
        ax=ax
    )
    
    ax.axhline(y=50, color='gray', linestyle='--', linewidth=1)
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 100)
    ax.set_title(f"Poly{aa}", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("")

# 不要な余分なサブプロット枠を削除（20個未満の場合用）
for j in range(len(amino_acids), len(axes)):
    fig.delaxes(axes[j])

# === 全体装飾 ===
fig.suptitle("Repeat Length vs Avg pLDDT for Each PolyX Type", fontsize=18)
#fig.text(0.5, 0.04, "Repeat Length", ha='center', fontsize=14)
#fig.text(0.08, 0.5, "Avg pLDDT", va='center', rotation='vertical', fontsize=14)
plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])

# === 保存 ===
plt.savefig("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig6/Where_is_PolyX_in_protein/scatter_all_polyX_25.png", dpi=300)
