import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"

#-----読み込むファイルの指定-----
TECAN_data = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/NI" #テカンデータを保存しているディレクトリを指定
samples = ['Δ','A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N','P', 'Q', 'R','S','T','V','W','Y'] #サンプル名前を指定
dirs = ["0", "50", "150", "200", "300", "400", "500_1", "500_2"]# ディレクトリ名のリスト
dirs_2 = ["0", "50", "150", "200", "300", "400", "500"]# ディレクトリ名のリスト
colors = {"0": "#AFEFF9",
          "50": "#71E7FA",
          "150": "#009BB3",
          "200": "#C0B1C9",
          "300":    "#D091FA",
          "400":    "#B241FD",
          "500":     "#6600AA"}
#-----データを整形する----------
dfs = {}
for d in dirs:
    file_path = os.path.join(TECAN_data, d, "maindata.xlsx")
    df = pd.read_excel(file_path, sheet_name="Data")
    if "Unnamed: 0" in df.columns:
        df.rename(columns={"Unnamed: 0": "sample"}, inplace=True)
    df["aTc (nM)"] = d
    dfs[d] = df
# 全てのデータを結合する
all_df = pd.concat(dfs, ignore_index=True)
all_df["aTc (nM)"] = all_df["aTc (nM)"].replace({"500_1": "500", "500_2": "500"})
#-----正規化を行う------------
# Vector（Φ）のみを抽出してMGRの平均を計算（Φ_1～Φ_4）
phi_rows = all_df[(all_df['aTc (nM)'] == '500') & (all_df['sample'].str.startswith('Δ_'))]
MGR100 = phi_rows['MGR'].mean()

# EGFP（Δ）のみを抽出してMFIの平均を計算（Δ_1～Δ_4）
delta_rows = all_df[(all_df['aTc (nM)'] == '500') & (all_df['sample'].str.startswith('Δ_'))]
MFI100 = delta_rows['MFI'].mean()

all_df['MGR_pct'] = all_df['MGR'] / MGR100 * 100 #　MGRを正規化
all_df['MFI_pct']  = all_df['MFI']  / MFI100  * 100 #　MFIを正規化
#all_df['MFI_pct'] = all_df['MFI_pct'].clip(lower=0.1) #　MFIが0.1以下のサンプルは0.1にする
all_df['NI_pct*pct']  = all_df['MGR_pct'] * all_df['MFI_pct'] #　NIを正規化

#-----グラフ化する----------------------
# サブプロットのサイズを調整（例：5行 × 4列）(8,10),（例：6行 × 4列）(8,12)
n = len(samples)
fig, axes = plt.subplots(5, 4, figsize=(8,10), sharex=True, sharey=True)
axes = axes.flatten()
max_ni_list = []

for ax, samp in zip(axes, samples):
    sub = all_df[all_df["sample"].str.startswith(samp + "_")]

    #最終的なNIを示す
    grouped = sub.groupby('aTc (nM)', as_index=False).mean(numeric_only=True)
    init_row = grouped.loc[grouped['NI_pct*pct'].idxmax()]# 濃度ごとに平均値を取り、NIの平均が最大となる濃度を採用
    #init_row = grouped.loc[grouped['MFI_pct'].idxmax()]# 濃度ごとに平均値を取り、MFIの平均が最大となる濃度を採用
    #init_row = sub.loc[sub['MFI_pct'].idxmax()]#　n=1でMFIが最大となる濃度を採用
    md = init_row["aTc (nM)"]
    max_row = sub[sub["aTc (nM)"] == md]

    max_ni = max_row['NI_pct*pct'].mean()
    max_ni_std = max_row['NI_pct*pct'].std()
    max_mfi = max_row['MFI_pct'].mean()
    max_mgr = max_row['MGR_pct'].mean()
    for d in dirs_2:
        dsub = sub[sub["aTc (nM)"] == d]
        ax.scatter(dsub["MFI_pct"], dsub["MGR_pct"], label=d, color=colors[d], s=20, alpha=0.8)
    # NIを表す面積を作図する（最大点に垂直・水平線を引く）
    ax.vlines(x=max_mfi, ymin=0, ymax=max_mgr, colors='gray', linestyles='--', linewidth=0.5)
    ax.hlines(y=max_mgr, xmin=0, xmax=max_mfi, colors='gray', linestyles='--', linewidth=0.5)
    ax.set_xscale('log')
    ax.set_xlim(0.05, 1000) 
    # タイトルにサンプル名とNI値
    ax.set_title(f"{samp} : {int(max_ni)}")
    ax.grid(False)
    # Excel用リストにNIデータを追加
    max_ni_list.append({"sample": samp, "aTc (nM)": md, "NI": max_ni, "SD": max_ni_std})

# 余分なサブプロットを消す
for ax in axes[n:]:
    fig.delaxes(ax)

fig.tight_layout()
plt.savefig(f"{TECAN_data}/NI.png", dpi=350)

ni_list = []
for samp in samples:
    row = {"sample": samp}
    for d in dirs_2:
        sub = all_df[all_df["sample"].str.startswith(samp + "_")]
        d_row = sub[sub["aTc (nM)"] == d]
        ni = d_row['NI'].mean()
        ni_std = d_row['NI'].std()
        row[d] = ni
        #row[f"{d}_std"] = ni_std
    ni_list.append(row)


#結果をExcelで出力
max_ni_df = pd.DataFrame(max_ni_list)
all_ni_df = pd.DataFrame(ni_list)
out_xlsx = os.path.join(TECAN_data, "NI_.xlsx")
with pd.ExcelWriter(out_xlsx) as writer:
    all_df.to_excel(writer, index=False, sheet_name="All_Data")# 全データ（正規化前＋正規化後すべて）
    all_ni_df.to_excel(writer, index=False, sheet_name="All_NI_Data")
    max_ni_df.to_excel(writer, index=False, sheet_name="Max_NI_pct") # 各サンプルの図中NI_pct

#-----グラフ化する----------------------
def create_NIplot(All_df, NI_df):
    NI_df = NI_df.sort_values(by="NI", ascending=False).reset_index(drop=True) #大きい順にソート
    fig, ax = plt.subplots(figsize=(8, 3))# プロットの準備

   # NI 平均値と標準偏差の棒グラフ
    bars = ax.bar(NI_df['sample'], NI_df['NI'], color=np.where(NI_df['sample'].isin(["Δ"]), "#7F7F7F", "#D0B9A4"), width=0.6, linewidth=0.5, edgecolor='black', alpha=0.7)
    ax.errorbar(NI_df['sample'], NI_df['NI'], yerr=NI_df['SD'], fmt='none', ecolor='black', elinewidth=0.5, capsize=0)

    # 生データのドットプロット
    for AA in samples:
        for d in dirs_2:
            available = All_df[All_df["sample"].str.startswith(f'{AA}_')]
            available = available[available["aTc (nM)"] == d]
            x_positions = [AA] * len(available)
            y_values_NI = available["NI_pct*pct"].values
            ax.scatter(x_positions, y_values_NI, color=colors[d], s=10, edgecolor='black', linewidths=0.1, zorder=3, alpha=0.5)
  
    # NIの範囲設定
    ax.set_ylim(1, 100000)
    ax.set_yscale('log')

    # プロットを保存
    plt.savefig(f"{TECAN_data}/NI_plot.png", dpi=350)

create_NIplot(all_df, max_ni_df)