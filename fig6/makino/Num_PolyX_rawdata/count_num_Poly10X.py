import re
from collections import Counter
import glob
import os
import pandas as pd

def count_genes_by_sequence_code(file_path):
    # ファイルの内容を読み込む
    with open(file_path, 'r') as file:
        content = file.read()

    # 正規表現で遺伝子IDとシーケンスコードを抽出
    gene_pattern = re.compile(r'^>([\w.-]+)\t([A-Z])', re.MULTILINE)
    matches = gene_pattern.findall(content)
    
    # 各シーケンスコードの最初の文字をカウント
    sequence_counts = Counter(sequence_code[1] for sequence_code in matches)
    return sequence_counts


# フォルダパスを指定
folder_path = '/Users/muraseyukihiro/Desktop/PolyX/Dry/makino/Num_PolyX_rawdata'

order = ['D', 'Q', 'N', 'S', 'E', 'R', 'H', 'P', 'K', 'A', 'V', 'T','G', 'L', 'C','F','I','Y','M','W']

# フォルダ内のすべての .fasta ファイルを取得して読み込む
data_list = []
for fasta_file in glob.glob(f"{folder_path}/*.fasta"):
    gene_counts = count_genes_by_sequence_code(fasta_file)
    gene_counts = {key: gene_counts[key] for key in order}
    df_counts = pd.DataFrame(list(gene_counts.items()), columns=['Key', os.path.basename(fasta_file).split('.')[0]])  # ファイル名をインデックスとして設定
    data_list.append(df_counts)

df = pd.concat(data_list, axis=1)

df.to_csv('/Users/muraseyukihiro/Desktop/PolyX/Dry/makino/Num_PolyX_rawdata/Num_PolyX_re.csv', index=False)



