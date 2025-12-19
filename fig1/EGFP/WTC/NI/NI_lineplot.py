import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smm
import os
plt.rcParams["font.family"] = "Arial"

#外れデータはサンプルシートでEmptyとして指定する

# ========================= ユーザー設定（入出力・描画） =========================
# --- 入力ファイル ---
TECAN_file_name_1 = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/0_500/20250623_WTCpro_EGFP_PolyX_aTc4ponts.xlsx" #TECANのraw_dataを指定
sample_sheet_name_1 = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/0_500/Sample_Sheet.xlsx" #サンプルシートを指定

TECAN_file_name_2 = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/200_500/20250626_WTCpro_EGFP_PolyX_aTc4points_re.xlsx" #TECANのraw_dataを指定
sample_sheet_name_2 = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/200_500/Sample_Sheet.xlsx" #サンプルシートを指定

# --- 出力先 ---
path = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig1/EGFP/WTC/NI"#保存先のパス

# --- サンプルの表示順（凡例・プロット順など）---
sample_name = ['Δ','A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N','P', 'Q', 'R','S','T','V','W','Y'] #サンプル名前を指定
MDs = ["_0_", "_50_", "_150_","_200_", "_300_", "_400_", "_500_" ]

# --- ラインプロット描画設定 ---
FIN_TIME_HOUR = 72        # x軸上限 [hour]
GROWTH_YLIM = (0, 1.5)    # 増殖曲線のylim
GROWTH_YTICKS = [0.5, 1.0, 1.5]
FI_YLIM = (0, 10000)     # GFP曲線のylim
FI_YTICKS = [4000, 8000]
PERIOD_MIN = 10

colors = {
    "_0_": "#AFEFF9",
    "_50_": "#71E7FA",
    "_150_": "#009BB3",
    "_200_": "#C0B1C9",
    "_300_": "#D091FA",
    "_400_": "#B241FD",
    "_500_": "#6600AA"
}
#-----------------------------

# ========================= データ読み込み・整形 =========================
TECAN_df_1 = pd.read_excel(TECAN_file_name_1, sheet_name=0, header=None)
sample_sheet_df_1 = pd.read_excel(sample_sheet_name_1, sheet_name=0)
well_to_sample_1 = dict(zip(sample_sheet_df_1['Well'], sample_sheet_df_1['Sample_number']))

TECAN_df_2 = pd.read_excel(TECAN_file_name_2, sheet_name=0, header=None)
sample_sheet_df_2 = pd.read_excel(sample_sheet_name_2, sheet_name=0)
well_to_sample_2 = dict(zip(sample_sheet_df_2['Well'], sample_sheet_df_2['Sample_number']))

# --- TECANのデータから各項目の範囲を読み込み、Sample_number列と結合する ---
def tecan_data(TECAN_df, well_to_sample, label):
    start_index = TECAN_df[TECAN_df[0].astype(str).str.contains(label)].index[0] + 4 #何行目から始まるのか
    end_index = TECAN_df.iloc[start_index:, 0].isna().idxmax()#何行目に終わるのか

    df = TECAN_df.iloc[start_index:end_index].copy()#範囲の指定
    df.dropna(how="any", axis=1, inplace=True)
    df.loc[:, 0] = df[0].map(well_to_sample)
    df.set_index(0, inplace=True)#Sample_numberを行名に指定
    if df.shape[1] % 2 != 0:  # 列数が奇数かどうかを確認
        df.drop(df.columns[-1], axis=1, inplace=True)
        print("サイクル数を偶数にする必要があるため、最後の列を削除しました。")
    df = df.apply(pd.to_numeric)#念の為数値データに直す
    return df

# --- 自家蛍光（Φ）の値を引く ---
def subtract_phi(fi_df: pd.DataFrame) -> pd.DataFrame:
    phi_samples = [s for s in ["Φ_1", "Φ_2", "Φ_3", "Φ_4"] if s in fi_df.index]
    if len(phi_samples) == 0:
        print("Φ サンプルが見つからないため、fiの差し引きをスキップします。")
        return fi_df
    phi_mean = fi_df.loc[phi_samples].mean(axis=0)  # 各時刻の平均
    fi_df = fi_df.subtract(phi_mean, axis=1).clip(lower=0)
    return fi_df

# ========================= 可視化（ラインプロット） =========================
def plot_filled_line(ax, AA, data, time, color, label, ylim=None, yticks=None):
    alpha=0.2
    linewidth=4
    mean = data.mean(axis=0)
    #std = data.std(axis=0)
    
    ax.plot(time, mean, color=color, linewidth=linewidth)
    #ax.fill_between(time, mean - std, mean + std, alpha=alpha, color=color, label=label)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_title(AA, fontsize="15")

def create_lineplot(df, ylim, yticks, lineplot_file):
    fig = plt.figure(dpi=100, figsize=(20, 10))
    time = np.arange(len(df.T)) / (60/PERIOD_MIN)  # interval time : 10min or 30min

    for n, AA in zip(range(len(sample_name)), sample_name):
        ax = fig.add_subplot(4, 6, n + 1)

        for md in MDs:
            data_md = df[df.index.str.startswith(f'{AA}{md}')]
            #data_md = df.loc[available]
            plot_filled_line(ax, AA, data_md, time, color=colors.get(md, "#23375B"), label="_nolegend_", ylim=ylim, yticks=yticks)

        # x軸の表示範囲
        ax.set_xlim(0, FIN_TIME_HOUR)
        ax.set_xticks([0, FIN_TIME_HOUR / 2, FIN_TIME_HOUR])

        # グラフのスタイル設定
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_linewidth(2)
        ax.spines["right"].set_linewidth(2)
        ax.spines["bottom"].set_linewidth(2)
        ax.tick_params(width=1, labelsize=12)

    plt.subplots_adjust(wspace=0.6, hspace=0.45)
    plt.savefig(lineplot_file, dpi=350)

#=================================================================================================================================


Growth_df_full_1 = tecan_data(TECAN_df_1, well_to_sample_1 ,'Label1')
Growth_df_full_2 = tecan_data(TECAN_df_2, well_to_sample_2 ,'Label1')
FI_df_full_1 = tecan_data(TECAN_df_1, well_to_sample_1 ,'Label2')
FI_df_full_2 = tecan_data(TECAN_df_2, well_to_sample_2 ,'Label2')

Growth_df_full = pd.concat([Growth_df_full_1, Growth_df_full_2]).dropna(axis=1)
FI_df_full = pd.concat([FI_df_full_1, FI_df_full_2]).dropna(axis=1)
print(Growth_df_full)
#------------------------------------
create_lineplot(Growth_df_full, GROWTH_YLIM, GROWTH_YTICKS, lineplot_file = f'{path}/growth.png')
create_lineplot(FI_df_full,FI_YLIM, FI_YTICKS, lineplot_file = f'{path}/fluorescence.png')