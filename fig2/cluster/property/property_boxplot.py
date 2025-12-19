import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==== 入力 ====
df_o = pd.read_excel('/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/AAindex1_property_18.xlsx')

x_order = ['Total', '1', '2', '3']

group_mapping = {
    'Total': ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N','P', 'Q', 'R','S','T','V','W','Y'],
    '1': ['E', 'S', 'N', 'Q'],
    '2': ['C', 'P', 'K', 'R', 'M', 'V', 'F', 'I','Y', 'L', 'W'],
    '3': ['G', 'T', 'H', 'A', 'D'],
}

# プロット対象
features = df_o.drop(columns=["AminoAcid"]).columns.tolist()

# グループ展開
rows = []
for _, row in df_o.iterrows():
    aa = row['AminoAcid']
    for g, members in group_mapping.items():
        if aa in members:
            rr = row.copy()
            rr['Group'] = g
            rows.append(rr)
df = pd.DataFrame(rows)

# 描画ループ
for feature in features:
    fig, ax = plt.subplots(figsize=(4, 3))
    
    sns.boxplot(
        data=df, x='Group', y=feature, order=x_order, width=0.3,
        showcaps=True, boxprops={'facecolor':'none'}, whiskerprops={'linewidth':1.5},
        medianprops={'color':'black'}, ax=ax, showfliers=False
    )
    sns.swarmplot(
        data=df, x='Group', y=feature, order=x_order,
        hue='AminoAcid', color="black", dodge=False, size=6, alpha=0.5, ax=ax
    )

    for i, row in df.iterrows():
        x = x_order.index(row['Group'])
        y = row[feature]
        ax.text(
            x, y, row['AminoAcid'],
            fontsize=6, ha='center', va='bottom', alpha=0.8
        )

    ax.legend_.remove()
    ax.set_title(feature, fontsize=5)
    ax.tick_params(labelsize=15)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.tight_layout()
    outpath = f"/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/cluster/property/{feature}.png"
    plt.savefig(outpath, dpi=350)
    plt.close(fig)