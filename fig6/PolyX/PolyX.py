import pandas as pd
import re

# CSVファイルを読み込む
file_path = '/Users/muraseyukihiro/Desktop/Sc_ORFtrans.csv'#SGDで出力したcsvファイルを指定(Systematic Name、Sequence列を入れておく)
data = pd.read_csv(file_path)
data = data[data['Qualifier'] != "Dubious"].reset_index(drop=True)#Dubiousを除く

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

# 各シーケンスに対して最大連続回数を計算し、新しいデータフレームを作成
max_counts_list = data['Sequence'].apply(PolyX_count)
max_counts_df = pd.DataFrame(max_counts_list.tolist())

#PolyXmaxとNum_Poly10Xをまとめたデータを作る
PolyXmax = max_counts_df.max()
Num_Poly10X = (max_counts_df >= 10).sum()
results_df = pd.DataFrame({'PolyXmax': PolyXmax,
                        'Num_Poly10X': Num_Poly10X})

max_counts_df.insert(0, 'Gene', data['Systematic name'])
max_counts_df.insert(1, 'Sequence', data['Sequence'])

# 結果をExcelファイルに保存する
output_file_path = '/Users/muraseyukihiro/Desktop/a.xlsx'
with pd.ExcelWriter(output_file_path, engine='openpyxl') as writer:
    max_counts_df.to_excel(writer, sheet_name='Max Counts', index=False)
    results_df.to_excel(writer, sheet_name='result')