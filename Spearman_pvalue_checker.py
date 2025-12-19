import pandas as pd
from scipy.stats import pearsonr, spearmanr

data = pd.read_excel('data.xlsx', sheet_name="raw", index_col=0)
print(data)

for i in data.columns:
    for j in data.columns:
        #ピアソン相関を求める
        p_correlation, p_pvalue = pearsonr(data[i], data[j])

        #スピアマンの順位相関を求める
        s_correlation, s_pvalue = spearmanr(data[i], data[j])

        print(f'{i} vs {j}')

        print(f"順位相関:{s_correlation}")
        print(f"順位相関:p={s_pvalue}")
        print(f"ピアソン相関:{p_correlation}")
        print(f"ピアソン相関:p={p_pvalue}")

        print("---------------")