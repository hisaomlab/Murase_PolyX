import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu

strains = ["BY4741", "cdc24", "mac1", "mmr1", "pre7", "psk1", "rpl18b", "rpl19a"]
plasmids = ["phi", "delta", "PolyE", "mox"]
directory_path = f'/Users/muraseyukihiro/Desktop/PolyX/Microscope/elongation_PolyE/all_data'  # CSVファイルが保存されているディレクトリ
column_name = ['AreaShape_MajorAxisLength', 'AreaShape_MinorAxisLength']    # 抽出したい列名
output_file = '/Users/muraseyukihiro/Desktop/PolyX/Microscope/elongation_PolyE/all_data/result.csv'  # 保存するファイル名

# ディレクトリ内のCSVファイルをすべて読み込み、特定の列を抽出してまとめる
all_data = []

# ディレクトリ内のすべてのCSVファイルを取得
for strain in strains:
    strain_path = os.path.join(directory_path, strain)
    for plasmid in plasmids:
        plasmid_path = os.path.join(strain_path, plasmid)
        for file in os.listdir(plasmid_path):
            if file == 'MyExpt_Segment.csv':
                file_path = os.path.join(plasmid_path, file)
                
                # CSVファイルを読み込み、特定の列を抽出
                df = pd.read_csv(file_path)
                extracted_data = df[column_name].copy()
                extracted_data['Major/Minor'] = (extracted_data['AreaShape_MajorAxisLength'] / extracted_data['AreaShape_MinorAxisLength'])
                extracted_data['strain'] = strain  # サンプル名を記録
                extracted_data['plasmid'] = plasmid  # サンプル名を記録
                all_data.append(extracted_data)

# すべてのデータを結合
combined_data = pd.concat(all_data, ignore_index=True)
# 結果をCSVに出力
combined_data.to_csv(output_file, index=False)

sample_colors = {'delta': 'gray', 'PolyE': '#73bf4d', 'phi': 'white', 'mox': 'yellow'} 

# strainごとにバイオリンプロットを作成
for strain in strains:
    strain_data = combined_data[combined_data['strain'] == strain]

    plt.figure(figsize=(10, 5))
    # 枠線を太くする
    ax = plt.gca()  # 現在の軸を取得
    for spine in ax.spines.values():
        spine.set_linewidth(2)  # 枠線の太さを設定
    sns.violinplot(x="plasmid", y="Major/Minor", data=strain_data, scale="count", palette=sample_colors)
    # 個別に輪郭線を調整
    for violin in plt.gca().collections:
        violin.set_edgecolor("black")  # 輪郭線を赤に変更
        violin.set_linewidth(1.5)    # 太さを調整
    #sns.swarmplot(x="plasmid", y="Major/Minor", data=strain_data, size=1.5, palette=sample_colors, edgecolor='black', linewidth=0.2)
    plt.ylim(1, 4)
    plt.axhline(y=1.5, color='black', linestyle='--', linewidth=1)  # 赤の点線、太さ2
    #plt.yscale('log')

    # プラズミドの組み合わせを設定
    plasmid_combinations = [("phi","delta"), ("phi", "PolyE"), ("phi", "mox")]
    # p値を計算しボンフェローニ補正を適用
    p_vals = []
    for plasmid1, plasmid2 in plasmid_combinations:
        #t検定
        #_, p_val = ttest_ind(strain_data[strain_data["plasmid"] == plasmid1]["Major/Minor"], strain_data[strain_data["plasmid"] == plasmid2]["Major/Minor"])
        #p_vals.append(p_val)
        #U検定
            _, p_val= mannwhitneyu(strain_data[strain_data["plasmid"] == plasmid1]["Major/Minor"], strain_data[strain_data["plasmid"] == plasmid2]["Major/Minor"])
            p_vals.append(p_val)
    # ボンフェローニ補正
    adjusted_p_values = [p * len(plasmid_combinations) for p in p_vals]

    # プロットのラベルを設定
    plt.title(f"Violin Plot of Major/Minor in {strain}", fontsize=16)
    plt.xlabel("plasmid", fontsize=12)
    plt.ylabel("Major/Minor", fontsize=12)
    plt.xticks(rotation=45, fontsize=10)

    # アノテーションを追加
    y_max = strain_data["Major/Minor"].max() * 1.2  # y軸の最大値に余裕を持たせる
    for i, (plasmid1, plasmid2) in enumerate(plasmid_combinations):
        if adjusted_p_values[i] < 0.05:
            pos1 = strain_data["plasmid"].unique().tolist().index(plasmid1)
            pos2 = strain_data["plasmid"].unique().tolist().index(plasmid2)
            #plt.plot([pos1, pos2], [y_max, y_max], color="black", lw=1.5)  # 線を引く
            #plt.text((pos1 + pos2) / 2, y_max * 1.05, f"p={adjusted_p_values[i]:.2e}", ha="center", va="bottom", fontsize=10)
            y_max *= 1.2  # 次のアノテーションの位置を少し上にする

    # プロットを表示
    plt.tight_layout()
    plt.savefig(f'/Users/muraseyukihiro/Desktop/PolyX/Microscope/elongation_PolyE/all_data/{strain}.pdf')
    # エクセルファイルに保存
    strain_data.to_excel(f'/Users/muraseyukihiro/Desktop/PolyX/Microscope/elongation_PolyE/all_data/data_{strain}.xlsx', index=False)

