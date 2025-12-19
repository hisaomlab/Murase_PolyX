import pandas as pd
from scipy.stats import pearsonr, spearmanr

data = pd.read_excel('/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/data.xlsx')
print(data)

a = ["S.cerevisiae_TDH3", "S.cerevisiae_WTC", "E.coli", "COS-7", "S.cerevisiae_U"]

for aa in a:
    for ab in a:
        #max_dataを使ってピアソン相関を求める
        p_correlation, p_pvalue = pearsonr(data[aa], data[ab])

        #max_dataを使ってスピアマンの順位相関を求める
        s_correlation, s_pvalue = spearmanr(data[aa], data[ab])

        print(f'{aa} vs {ab}')

        print(f"ピアソン相関:{p_correlation}")
        print(f"ピアソン相関:p={p_pvalue}")
        print(f"順位相関:{s_correlation}")
        print(f"順位相関:p={s_pvalue}")
        print("---------------")