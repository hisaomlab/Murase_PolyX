import pandas as pd
import re
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
input_dir = BASE_DIR / 'ClusterFasta'

if not input_dir.is_dir():
    raise FileNotFoundError(f'Input directory not found: {input_dir}')

records = []
fasta_suffixes = {'.fa', '.fasta'}

def iter_fasta_records(fasta_path: Path):
    header_parts = None
    seq_lines = []
    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header_parts is not None:
                    yield header_parts[0], ''.join(seq_lines)
                header_parts = line[1:].split()
                seq_lines = []
            else:
                seq_lines.append(line)
    if header_parts is not None:
        yield header_parts[0], ''.join(seq_lines)

fasta_files = sorted(
    path for path in input_dir.rglob('*')
    if path.is_file() and path.suffix.lower() in fasta_suffixes
)

for fasta_file in fasta_files:
    for name, sequence in iter_fasta_records(fasta_file):
        records.append({
            'Systematic name': name,
            'Sequence': sequence
        })

if not records:
    raise ValueError(f'No FASTA records found under {input_dir}')

data = pd.DataFrame(records)

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

# 詳細結果はCSV、集計結果はExcelで保存する
max_counts_csv_path = BASE_DIR / 'max_counts.csv'
max_counts_df.to_csv(max_counts_csv_path, index=False)

output_file_path = BASE_DIR / 'sammary.xlsx'
results_df.to_excel(output_file_path, sheet_name='result')
