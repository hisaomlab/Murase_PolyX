import pandas as pd
from pandas.api.types import is_numeric_dtype
from scipy.stats import spearmanr

# Excelファイル読み込み
data = pd.read_excel(
    '/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/AAindex1_property_18.xlsx',
    sheet_name='AAindex1_property_18'
)

NI_lists = ["EGFP", "moxGFP", "mNeonGreen", "Gamillus", "mScarlet-I", "mCherry"]

# 計算対象プロパティ列（AminoAcid と NI 列は除外、かつ数値列のみ）
exclude_cols = set(['AminoAcid'] + NI_lists)
property_cols = [c for c in data.columns if c not in exclude_cols and is_numeric_dtype(data[c])]
results = []

for NI in NI_lists:
    if NI not in data.columns:
        continue  # 列が存在しない場合はスキップ

    ni_series = pd.to_numeric(data[NI], errors='coerce')

    for prop in property_cols:
        prop_series = pd.to_numeric(data[prop], errors='coerce')
        valid = pd.concat([ni_series, prop_series], axis=1).dropna()

        if len(valid) < 2:
            s_corr = float('nan')
        else:
            s_corr, _ = spearmanr(valid[NI], valid[prop])

        results.append({
            "Property": prop,
            NI: round(s_corr, 4) if pd.notna(s_corr) else s_corr
        })

# Property ごとにまとめる（横に NI 列が並ぶ形式）
results_df = pd.DataFrame(results).groupby("Property", as_index=False).first()

# Excel出力（Spearman のみ）
output_path = '/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/AAindex/18/NI_correlation_results.xlsx'
results_df.to_excel(output_path, index=False, sheet_name='Spearman')

print(f"✅ Spearman相関結果をExcelに出力しました: {output_path}")