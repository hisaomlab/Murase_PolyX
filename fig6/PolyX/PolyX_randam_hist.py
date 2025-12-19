import pandas as pd
import re
import random
import matplotlib.pyplot as plt
import numpy as np

# CSVファイルを読み込む
file_path = '/Users/muraseyukihiro/Desktop/PolyX/ORF2024109.csv'
data = pd.read_csv(file_path)
data = data[data['Qualifier'] != "Dubious"].reset_index(drop=True)#Dubiousを除く

# アミノ酸配列をランダムに入れ替える関数
def shuffle_sequence(sequence):
    sequence = sequence.replace('*', '')# アスタリスクを取り除く
    composition = list(sequence)
    random.shuffle(composition)
    return ''.join(composition) + '*'

# 各アミノ酸の最大連続回数を計算する関数
def PolyX_count(sequence):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    max_counts = {}
    for aa in amino_acids:
        try:
            max_counts[aa] = max(len(match) for match in re.findall(re.escape(aa) + "+", sequence))
        except ValueError:
            max_counts[aa] = 0
    return max_counts


shuffle_PolyXmax_list = []
shuffle_Num_Poly10X_list = []

for i in range(10000):
    SEED = i
    random.seed(SEED)

    data['Shuffled_Sequence'] = data['Sequence'].apply(shuffle_sequence)# シーケンスをランダムに入れ替える

    max_counts_list = data['Shuffled_Sequence'].apply(PolyX_count)# 各シーケンスに対して最大連続回数を計算し、新しいデータフレームを作成
    max_counts_df = pd.DataFrame(max_counts_list.tolist())

    #PolyXmaxとNum_Poly10Xをまとめたデータを作る
    PolyXmax = max_counts_df.max()
    Num_Poly10X = (max_counts_df >= 10).sum()
    shuffle_PolyXmax_list.append(PolyXmax)
    shuffle_Num_Poly10X_list.append(Num_Poly10X)
    print(i)

# リストをデータフレームに変換
shuffle_PolyXmax_df = pd.concat(shuffle_PolyXmax_list, axis=1)
shuffle_Num_Poly10X_df = pd.concat(shuffle_Num_Poly10X_list, axis=1)

print(shuffle_PolyXmax_df)

mean1 = shuffle_PolyXmax_df.mean(axis=1)
max1 = shuffle_PolyXmax_df.max(axis=1)
min1 = shuffle_PolyXmax_df.min(axis=1)
std1 = shuffle_PolyXmax_df.std(axis=1)
PolyXmax_summary_df = pd.DataFrame({'Min': min1,
                           'Mean': mean1,
                           'Max': max1,
                           'SD': std1})

mean2 = shuffle_Num_Poly10X_df.mean(axis=1)
max2 = shuffle_Num_Poly10X_df.max(axis=1)
min2 = shuffle_Num_Poly10X_df.min(axis=1)
std2 = shuffle_Num_Poly10X_df.std(axis=1)
Num_Poly10X_summary_df = pd.DataFrame({'Min': min2,
                           'Mean': mean2,
                           'Max': max2,
                           'SD': std2})

#結果をExcelファイルに保存する
output_file_path = '/Users/muraseyukihiro/Desktop/PolyX/Dry/PolyX/PolyX_shuffle_analysis_10000.xlsx'
with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
    shuffle_PolyXmax_df.to_excel(writer, sheet_name='shuffle_PolyXmax')
    shuffle_Num_Poly10X_df.to_excel(writer, sheet_name='shuffle_Num_Poly10X')
    PolyXmax_summary_df.to_excel(writer, sheet_name='PolyXmax_result')
    Num_Poly10X_summary_df.to_excel(writer, sheet_name='Num_Poly10X_result')

# サブプロットの行数と列数を設定
num_plots = 20  # 表示したいヒストグラムの数
rows = 5  # 行数
cols = 4  # 列数

# 一枚の図を作成
fig, axes = plt.subplots(rows, cols, figsize=(15, 10))
axes = axes.flatten()  # 2次元配列を1次元に変換

# ヒストグラムを描く
for i in range(20):
    row_data = shuffle_PolyXmax_df.iloc[i]
    axes[i].hist(row_data, bins=np.arange(0.5, 51.5, 1), alpha=0.5)
    axes[i].set_xlim(1, 50)
    axes[i].set_xticks(np.arange(0, 51, 5))
    axes[i].set_title(shuffle_PolyXmax_df.index[i])
    axes[i].set_xlabel('PolyXmax')
    axes[i].set_ylabel('Frequency')

# 不要な空白のサブプロットを削除
for j in range(num_plots, len(axes)):
    fig.delaxes(axes[j])

# レイアウトを調整
plt.tight_layout()

# 図をファイルとして保存
plt.savefig('/Users/muraseyukihiro/Desktop/PolyX/Dry/PolyX/histograms.png')

# グラフを表示
plt.show()

