import pandas as pd
from scipy.stats import ttest_ind

# 読み込むExcelの設定
EXCEL_PATH = '/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig5/PolyE_protein/p_val.xlsx'
SHEET_NAME = "Sheet1"
TARGET_GROUPS = {"EGFP", "PolyE"}
TARGET_COLUMNS = ["_All", "_Sol", "_Insol"]
EQUAL_VAR = False  # True: Studentのt検定, False: Welchのt検定

# データを整形
raw = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_NAME)
print(raw)

# t検定を実行しp値を出力
for frac in TARGET_COLUMNS:
    egfp = raw[f"EGFP{frac}"].dropna().to_numpy(dtype=float)
    polye = raw[f"PolyE{frac}"].dropna().to_numpy(dtype=float)
    if len(egfp) < 2 or len(polye) < 2:
        print(f"{frac}: サンプル数が不足しています")
        continue
    _, p_value = ttest_ind(egfp, polye, equal_var=EQUAL_VAR)
    print(f"{frac}: p-value = {p_value:.4g}")
