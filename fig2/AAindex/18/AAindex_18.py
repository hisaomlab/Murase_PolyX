import re
import csv
from collections import defaultdict
import math

file = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/aaindex1"
# 正しいアミノ酸順（1文字）
AA_order = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
]

# 抜き出したいAAindex ID
target_ids = {
    'BUNA790103', 'FINA910104', 'GEOR030103', 'GEOR030104', 'LEVM760103',
    'MITS020101', 'NADH010107', 'NAKH920107', 'PALJ810107', 'QIAN880138',
    'RICJ880104', 'RICJ880117', 'ROBB760107', 'TANS770102', 'TANS770108',
    'VASM830101', 'WERD780103', 'WOEC730101'
}

# アミノ酸ごとの特性値を記録する辞書
aa_properties = {aa: {} for aa in AA_order}

with open(file, "r") as f:
    content = f.read()

# エントリごとに分割
entries = content.strip().split('//')

for entry in entries:
    lines = entry.strip().split('\n')
    if not lines:
        continue

    # アクセッション番号の取得
    acc_line = next((line for line in lines if line.startswith("H ")), None)
    if not acc_line:
        continue
    acc_id = acc_line[2:].strip()
    if acc_id not in target_ids:
        continue  # 指定されたIDでない場合はスキップ

    # Description 行
    desc_line = next((line for line in lines if line.startswith("D ")), None)
    if not desc_line:
        continue
    description = desc_line[2:].strip()

    # "I" 行から値を取得
    try:
        i_index = next(i for i, line in enumerate(lines) if line.startswith("I"))
        values_lines = lines[i_index+1:i_index+4]
    except StopIteration:
        continue

    values = []
    for line in values_lines:
        for val in line.strip().split():
            if val.upper() == 'NA':
                values.append(float('nan'))  # または 'NA' として文字列で残すなら str('NA')
            else:
                try:
                    values.append(float(val))
                except ValueError:
                    values.append(float('nan'))

    if len(values) != 20:
        continue  # 無効なエントリをスキップ

    for aa, val in zip(AA_order, values):
        aa_properties[aa][description] = val

# すべての特性を収集して列順を決定
all_descriptions = sorted({desc for props in aa_properties.values() for desc in props})

# CSV 出力
with open("/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/AAindex1_property_18.csv", "w", newline='') as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["AminoAcid"] + all_descriptions)

    for aa in AA_order:
        row = [aa]
        for desc in all_descriptions:
            val = aa_properties[aa].get(desc, "")
            row.append("" if (val is None or (isinstance(val, float) and math.isnan(val))) else val)
        writer.writerow(row)

print("✅ 出力完了：aaindex1_wide.csv")