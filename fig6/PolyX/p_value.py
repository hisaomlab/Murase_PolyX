import pandas as pd
import numpy as np

# ====== 入力 ======
shuffle_path = '/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig6/PolyX/PolyX_shuffle_analysis_10000.xlsx'
real_path    = '/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig6/PolyX/PolyX_analysis.xlsx'
real_sheet, shuffle_sheet = "result", "shuffle_PolyXmax"

# ====== 読み込み ======
real_df = pd.read_excel(real_path, sheet_name=real_sheet).copy()
if 'AA' not in real_df.columns:
    real_df = real_df.rename(columns={'Unnamed: 0': 'AA'})
real_df = real_df[['AA','PolyXmax']].dropna()
real_df['AA'] = real_df['AA'].astype(str)

shuffle_df = pd.read_excel(shuffle_path, sheet_name=shuffle_sheet).copy()

# --- AA列を特定してindex化 ---
# 画像の通りなら先頭列がAA。名前が無くてもOK。
aa_col = shuffle_df.columns[0]
shuffle_df = shuffle_df.rename(columns={aa_col: 'AA'}).set_index('AA')

# --- すべての列を数値化（文字はNaNになる） ---
for c in shuffle_df.columns:
    shuffle_df[c] = pd.to_numeric(shuffle_df[c], errors='coerce')

# ====== 補助関数 ======
def right_tail_p(obs, null):
    return (np.sum(null >= obs) + 1) / (len(null) + 1)

def two_sided_p(obs, null):
    med = np.median(null)
    return (np.sum(np.abs(null - med) >= abs(obs - med)) + 1) / (len(null) + 1)

def bh_fdr(pvals):
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = q
    return np.clip(out, 0, 1)

# ====== 計算 ======
rows = []
for aa, obs in real_df[['AA','PolyXmax']].itertuples(index=False):
    # 対応する行を取得して数値ベクトル化（NaN除去）
    if aa in shuffle_df.index:
        null = shuffle_df.loc[aa].astype(float).dropna().to_numpy()
    else:
        # 念のため位置対応（AA名が一致しなかった場合）
        idx = real_df.index[real_df['AA'] == aa][0]
        null = shuffle_df.iloc[idx].astype(float).dropna().to_numpy()

    mean, sd = float(np.mean(null)), float(np.std(null, ddof=1))
    z = (obs - mean) / sd if sd > 0 else np.nan
    pct = (np.sum(null <= obs) / len(null)) * 100.0
    p_r = right_tail_p(obs, null)
    p_2 = two_sided_p(obs, null)

    rows.append([aa, obs, mean, np.median(null), sd, pct, z, p_r, p_2])

out = pd.DataFrame(rows, columns=[
    'AA','Observed','NullMean','NullMedian','NullSD',
    'Percentile(%)','Zscore','P_right(>=)','P_two_sided'
]).sort_values('Observed', ascending=False)

out['q_right(BH)'] = bh_fdr(out['P_right(>=)'].values)
out['q_two(BH)']   = bh_fdr(out['P_two_sided'].values)

out.to_excel('/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig6/PolyX/PolyX_shuffle_statistics.xlsx', 
             index=False)

# 表示を少し丸める
for c in ['NullMean','NullMedian','NullSD','Percentile(%)','Zscore','P_right(>=)','P_two_sided','q_right(BH)','q_two(BH)']:
    out[c] = out[c].astype(float).round(4)

print(out.to_string(index=False))
