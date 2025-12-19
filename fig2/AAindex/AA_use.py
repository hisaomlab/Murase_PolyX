import pandas as pd
from collections import Counter

# CSVファイルを読み込む
file_path = '/Users/muraseyukihiro/Desktop/PolyX/ORF2024109.csv'#SGDで出力したcsvファイルを指定(Systematic Name、Sequence列を入れておく)
data = pd.read_csv(file_path)
data = data[data['Qualifier'] != "Dubious"].reset_index(drop=True)#Dubiousを除く

# 全ての配列を1つの文字列に結合
all_sequences = ''.join(data['Sequence']).replace('*', '')

# 各アミノ酸の使用回数をカウント
amino_acid_counts = Counter(all_sequences)

# 結果を表示
print(amino_acid_counts)

# 結果をDataFrameに変換して表示
amino_acid_counts_df = pd.DataFrame(list(amino_acid_counts.items()), columns=['Amino Acid', 'UseCount'])
amino_acid_counts_df = amino_acid_counts_df.sort_values(by='Amino Acid', ascending=True).reset_index(drop=True)

# データフレームを表示
print(amino_acid_counts_df)
print(len(all_sequences))
amino_acid_counts_df.to_csv("/Users/muraseyukihiro/Desktop/PolyX/Dry/PolyX/amino_acid_counts.csv", index=False)