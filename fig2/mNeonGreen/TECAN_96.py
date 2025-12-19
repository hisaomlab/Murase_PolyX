"""
TECAN 解析スクリプト — 処理フロー
1) 入力データ（TECAN 生データ・サンプルシート）を読み込む
2) Growth_df（OD）・FI_df（蛍光）を作成し、Φを補正
3) ラインプロットで時系列の平均±SDを描画
4) MGR・MFI・NI を計算（NIはΔ基準で正規化）
5) 平均・SD・t検定（Bonferroni補正）を集計
6) Excel と PDF 図に出力
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smm
plt.rcParams["font.family"] = "Arial"

#外れデータはサンプルシートでEmptyとして指定する

# ========================= ユーザー設定（入出力・描画） =========================
# --- 入力ファイル ---
TECAN_file_name = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/mNeonGreen/20241217_mNeonGreen_PolyX.xlsx" #TECANのraw_dataを指定
sample_sheet_name = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/mNeonGreen/Sample_Sheet_Poly10X.xlsx" #サンプルシートを指定

# --- 出力先 ---
path = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/mNeonGreen"#保存先のパス
plate_visualization_file = f'{path}/platedata.xlsx' #最大蛍光値・MGRを保存するファイル
lineplot_file = f'{path}/lineplot.png' #lineplotを保存するファイル
mainplot_file = f'{path}/mainplot.png' #mainplotを保存するファイル
NIplot_file = f'{path}/NIplot.png'#NIplotを保存するファイル
maindata_file = f'{path}/maindata.xlsx'#maindataを保存するファイル

# --- サンプルの表示順（凡例・プロット順など）---
sample_name = ['Δ','A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N','P', 'Q', 'R','S','T','V','W','Y'] #サンプル名前を指定

# --- ラインプロット描画設定 ---
FIN_TIME_HOUR = 80        # x軸上限 [hour]
GROWTH_YLIM = (0, 1.5)    # 増殖曲線のylim
GROWTH_YTICKS = [0.5, 1.0, 1.5]
FI_YLIM = (0, 15000)     # GFP曲線のylim
FI_YTICKS = [5000, 10000, 15000]

# --- MGR 計算パラメータ ---
PERIOD_MIN = 10           # 測定間隔（min）
RANGE = 25                # 回帰レンジ（サイクル数）
STEP_SIZE = 10            # ウインドウの移動幅（大きいほど高速）
R2_THRESHOLD = 0.9        # R^2のしきい値
GR_START_CYCLE = 50       # GR計算開始位置（サイクル数）

# --- プロット描画設定 ---
FIGSIZE = (10, 3)          # 図のサイズ
MGR_YLIM = [0, 5]              # MGRのylim
MFI_YLIM = [0, 15000]           # MFIのylim
NI_YLIM = [0, 20000]     # NIのylim
TITLE = "mNeonGreen-PolyX (SC-LU)"      # 図のタイトル
#-----------------------------

# ========================= データ読み込み・整形 =========================
TECAN_df = pd.read_excel(TECAN_file_name, sheet_name=0, header=None)
sample_sheet_df = pd.read_excel(sample_sheet_name, sheet_name=0)
well_to_sample = dict(zip(sample_sheet_df['Well'], sample_sheet_df['Sample_number']))

# --- TECANのデータから各項目の範囲を読み込み、Sample_number列と結合する ---
def tecan_data(label):
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

# ---（必要時のみ使用）移動平均を取る---
def fourier_transform(Growth_df, FI_df):
    #~~~~パラメータ~~~~
    Start = 35 #周期検出の最初
    End = 192 #周期検出の最後
    Threshold = 0.1 #閾値を設定
    Max_Cycle = 30 #最大周期を設定
    #周期を計算し、移動平均を適用することでデータを滑らかにする-------------------------------------------------------
    cycle_df = Growth_df.T.loc[Start:End]
    all_period = []
    for column_name in cycle_df:
        f = cycle_df[column_name].values
        t = np.arange(len(f))

        assert t.size == f.size  # 時間軸の長さとデータの長さが同じであることを確認する
        assert np.unique(np.diff(t)).size == 1  # 時間間隔が全て一定であることを確認する

        T = (t[1] - t[0]) * t.size
        period = 1.0 / (np.arange(t.size / 2)[1:] / T)

        # パワースペクトル密度を計算
        f = f - np.average(f)         # 平均をゼロに。
        F = fftpack.fft(f)                          # 高速フーリエ変換
        Po = np.abs(F[1:(t.size // 2)]) ** 2 / T

        dominant_periods = [x for x in period[Po > Threshold] if x <= Max_Cycle]

        all_period.extend(dominant_periods)
    period = round(np.mean(all_period))
    print('周期は',period,'です。')

    Growth_ave_df =  Growth_df.rolling(period,axis=1,center=True,min_periods=1).mean()
    FI_ave_df =  FI_df.rolling(period,axis=1,center=True,min_periods=1).mean()

    return Growth_ave_df, FI_ave_df
    
# ========================= 可視化（ラインプロット） =========================
def plot_filled_line(ax, AA, data, time, color, label, ylim=None, yticks=None):
    alpha=0.2
    linewidth=4
    mean = data.mean(axis=0)
    std = data.std(axis=0)
    
    ax.plot(time, mean, color=color, linewidth=linewidth)
    ax.fill_between(time, mean - std, mean + std, alpha=alpha, color=color, label=label)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_title(AA, fontsize="15")

def create_lineplot(growth_df, fi_df):
    fig = plt.figure(dpi=100, figsize=(20, 10))
    time = np.arange(len(growth_df.T)) / (60/PERIOD_MIN)  # interval time : 10min or 30min

    for n, AA in zip(range(len(sample_name)), sample_name):
        ax = fig.add_subplot(4, 6, n + 1)
        available_samples = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in growth_df.index]#_1,_2,_3,_4の中であるものを取得する

        # 増殖曲線グラフ
        growth_data = growth_df.loc[available_samples]
        plot_filled_line(ax, AA, growth_data, time, color="#23375B", label="_nolegend_", ylim=GROWTH_YLIM, yticks=GROWTH_YTICKS)

        # 蛍光強度のグラフ
        ax2 = ax.twinx()
        fi_data = fi_df.loc[available_samples]
        plot_filled_line(ax2, AA, fi_data, time, color="#58C9A7", label="_nolegend_", ylim=FI_YLIM, yticks=FI_YTICKS)

        # x軸の表示範囲
        ax.set_xlim(0, FIN_TIME_HOUR)
        ax.set_xticks([0, FIN_TIME_HOUR / 2, FIN_TIME_HOUR])

        # グラフのスタイル設定
        for axis in [ax, ax2]:
            axis.spines["top"].set_visible(False)
            axis.spines["left"].set_linewidth(2)
            axis.spines["right"].set_linewidth(2)
            axis.spines["bottom"].set_linewidth(2)
            axis.tick_params(width=1, labelsize=12)

    plt.subplots_adjust(wspace=0.6, hspace=0.45)
    plt.savefig(lineplot_file, dpi=350)

#============================ MGR MFI NI計算関数 =====================================
# ---- MGR MFI NI算出 ----
def calculation(Growth_df: pd.DataFrame, FI_df: pd.DataFrame):
    #データをログ変換
    OD_log2 = np.log2(Growth_df).dropna()
    Max_cycle = len(OD_log2.columns)
    x = np.arange(RANGE)
    MGR = OD_log2[1]

    #WellごとにTimeウインドウを移動しながらGR（直線回帰の傾き=SlopeとR2を取得）、R2のThreshold以上のSlopeを返す（Threshold以下だと0.01を返す）。
    #最終的なGRは測定時間間隔（Period）で割り(min-1)、便宜的に1000をかけている。

    for index in MGR.index:
        Cycle = GR_START_CYCLE
        n = 1
        Slope  = np.array([0.1])
        R2 = np.array([1.0])
        OD_list = OD_log2.loc[index].values

        while Cycle < Max_cycle-RANGE:
            OD_list2 = OD_list[Cycle:Cycle+RANGE]
            Slope = np.insert(Slope,n,(np.polyfit(x,OD_list2,1)[0] * 1000 / PERIOD_MIN)) #1次関数で傾き（Slope）を計算。1000倍してPeriod(min)で割り、Slopeに順次代入
            R2 = np.insert(R2,n,(np.corrcoef(x,OD_list2)[0][1]**2)) #R2を計算。R2に順次代入
            Cycle += STEP_SIZE
            n += 1

        GR = np.array([x for x in Slope[R2 > R2_THRESHOLD]]) #R2がThresholdより高いGRのarrayを作成
        MGR[index] = GR.max()

    # ---- MFIの算出（0の場合は0.1にする） ----
    MFI = FI_df.max(axis=1).clip(lower=0.1)

    # ---- NIの算出（%MGR(/Δ), %MFI(/Δ) を計算 → NI = %MGR × %MFI) ----
    base_samples = [s for s in ["Δ_1", "Δ_2", "Δ_3", "Δ_4"] if s in MGR.index]
    if len(base_samples) == 0:
        raise ValueError("Δ のレプリケートが見つかりません（Δ_1〜Δ_4 のいずれかが必要）")

    mgr_base = MGR.loc[base_samples].mean()
    mfi_base = MFI.loc[base_samples].mean()
    # 0割回避（通常は起こらないが保険）
    if mgr_base == 0: mgr_base = 1e-9
    if mfi_base == 0: mfi_base = 1e-9

    MGR_pct = (MGR / mgr_base) * 100.0      # %MGR(/Δ)
    MFI_pct = (MFI / mfi_base) * 100.0      # %MFI(/Δ)
    NI = MGR_pct * MFI_pct                  # NI（単位は %^2）
    return MGR, MFI, NI, MGR_pct, MFI_pct

# ---- データの整形 ----
def create_maindata(MGR, MFI, NI):
    # 生データをまとめる
    df = pd.DataFrame({"MGR": MGR, "MFI": MFI, "NI": NI})
    mean_dict, sd_dict, p_raw_dict = {}, {}, {}

    # p値を出す際の基準のサンプルを取得
    base_samples = [s for s in ["Δ_1", "Δ_2", "Δ_3", "Δ_4"] if s in df.index]

    for AA in sample_name:
        available = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in df.index]#_1,_2,_3,_4の中であるものを取得する
        
        # 平均値
        mean_dict[AA] = df.loc[available].mean()
        
        # SD
        sd_dict[AA]   = df.loc[available].std()
        
        # t検定
        if base_samples:
            p_mgr = ttest_ind(df.loc[available, "MGR"], df.loc[base_samples, "MGR"], equal_var=True)[1]
            p_mfi = ttest_ind(df.loc[available, "MFI"], df.loc[base_samples, "MFI"], equal_var=True)[1]
            p_ni  = ttest_ind(df.loc[available, "NI"],  df.loc[base_samples, "NI"],  equal_var=True)[1]
        else:
            p_mgr = p_mfi = p_ni = np.nan

        p_raw_dict[AA] = {"MGR": p_mgr, "MFI": p_mfi, "NI": p_ni}
    
    df_mean = pd.DataFrame(mean_dict).T
    df_sd   = pd.DataFrame(sd_dict).T
    df_praw = pd.DataFrame(p_raw_dict).T 

    # ボンフェロー二補正
    if not df_praw.empty:
        p_mgr_adj = smm.multipletests(df_praw["MGR"].fillna(1.0), alpha=0.05, method='bonferroni')[1]
        p_mfi_adj = smm.multipletests(df_praw["MFI"].fillna(1.0), alpha=0.05, method='bonferroni')[1]
        p_ni_adj  = smm.multipletests(df_praw["NI"].fillna(1.0),  alpha=0.05, method='bonferroni')[1]
        df_p = pd.DataFrame({
            "MGR_bonferroni": p_mgr_adj,
            "MFI_bonferroni": p_mfi_adj,
            "NI_bonferroni":  p_ni_adj
        }, index=df_praw.index)
    else:
        df_p = pd.DataFrame(columns=["MGR_bonferroni", "MFI_bonferroni", "NI_bonferroni"])

    # Excelにデータフレームを書き込む
    with pd.ExcelWriter(maindata_file) as writer:
        df.to_excel(writer, sheet_name="Data")
        df_mean.to_excel(writer, sheet_name="Mean")
        df_sd.to_excel(writer, sheet_name="SD")
        df_p.to_excel(writer, sheet_name="P-Values")
    return {"Data": df, "Mean": df_mean, "SD": df_sd, "P-Values": df_p}

#============================ 可視化（プレート） =====================================
def plate_visualization(MGR, MFI, NI):
    #MGRの整形
    MGR = np.array(MGR.values)
    MGR = MGR.reshape([8,12])

    #MFIの整形
    MFI = np.array(MFI.values)
    MFI = MFI.reshape([8,12])

    #NIの整形
    NI = np.array(NI.values)
    NI = NI.reshape([8,12])

    column = [i for i in range(1, 13)]
    row = ["A", "B", "C", "D", "E", "F", "G", "H"]

    MGR = pd.DataFrame(MGR,columns=column, index=row)
    MFI = pd.DataFrame(MFI,columns=column, index=row)
    NI = pd.DataFrame(NI,columns=column, index=row)

    with pd.ExcelWriter(plate_visualization_file) as writer:
        MGR.to_excel(writer,sheet_name="MGR")
        MFI.to_excel(writer,sheet_name="MFI")
        NI.to_excel(writer,sheet_name="NI")

#============================ 可視化（MGR MFIプロット, NIプロット） =====================================
def create_mainplot(results):
    df      = results["Data"]
    df_mean = results["Mean"]
    df_sd   = results["SD"]
    df_p    = results["P-Values"]
    p_mgr_adj = df_p.get("MGR_bonferroni", pd.Series(1.0, index=df_mean.index)).reindex(df_mean.index)
    p_mfi_adj = df_p.get("MFI_bonferroni", pd.Series(1.0, index=df_mean.index)).reindex(df_mean.index)

    fig, ax = plt.subplots(figsize=FIGSIZE)# プロットの準備

   # MGR 平均値と標準偏差の棒グラフ
    bars = ax.bar(df_mean.index, df_mean['MGR'], color=np.where(df_mean.index.isin(["Φ"]), "gray", "#23375B"), label='Max growth rate', width=0.6, linewidth=0.5, edgecolor='black', zorder=1, alpha=0.7)
    ax.errorbar(df_mean.index, df_mean['MGR'], yerr=df_sd['MGR'], fmt='none', ecolor='black', elinewidth=0.5, capsize=0, zorder=2)

    # MFI 平均値と標準偏差のプロット
    ax2 = ax.twinx()
    ax2.scatter(df_mean.index, df_mean['MFI'], color='#58C9A7', marker='_', label='Expression level', s=150, edgecolor='black', linewidth=2.5, zorder=2, alpha=0.9)
    ax2.errorbar(df_mean.index, df_mean['MFI'], yerr=df_sd['MFI'], fmt='none', ecolor='black', elinewidth=0.5, capsize=0, zorder=1)

    # 生データのドットプロット
    for AA in sample_name:
        available = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in df.index]
        x_positions = [AA] * len(available)

        y_values_MGR = df.loc[available, "MGR"].values#MGRの生データ
        ax.scatter(x_positions, y_values_MGR, color='black', s=10, linewidths=0.1, zorder=3, alpha=0.5)

        y_values_MFI = df.loc[available, "MFI"].values#MFIの生データ
        ax2.scatter(x_positions, y_values_MFI, color='#58C9A7', s=10, linewidths=0.1, zorder=3, alpha=0.7)
  
    # MGRの範囲設定
    ax.set_ylim(MGR_YLIM)  # MGRの範囲で設定
    ax2.set_ylim(MFI_YLIM)  # MFIの範囲で設定

    # 有意差のあるサンプルについて*マークをつける処理(MGR)
    for i, (p_val, rect) in enumerate(zip(p_mgr_adj, bars)):
        if p_val < 0.05:  # 有意差のある場合
            ax.text(rect.get_x() + rect.get_width() / 2., rect.get_y() - 0.75, '*', ha='center', color='k', fontsize=10)
    
    # 有意差のあるサンプルについて*マークをつける処理(MFI)
    for i, (p_val, rect) in enumerate(zip(p_mfi_adj, bars)):
        if p_val < 0.05:  # 有意差のある場合
            ax.text(rect.get_x() + rect.get_width() / 2., rect.get_y() - 0.8, '*', ha='center', color='#58C9A7', fontsize=10)
  
    # ラベルとタイトルの設定
    ax.set_title(f'MGR and MFI: {TITLE}', fontsize=14)

    # プロットを表示
    #plt.show()
    plt.savefig(mainplot_file, dpi=350)

def create_NIplot(results):
    df      = results["Data"]
    df_mean = results["Mean"]
    df_sd   = results["SD"]
    df_p    = results["P-Values"]
    p_ni_adj = df_p.get("NI_bonferroni", pd.Series(1.0, index=df_mean.index)).reindex(df_mean.index)

    fig, ax = plt.subplots(figsize=FIGSIZE)# プロットの準備

   # NI 平均値と標準偏差の棒グラフ
    bars = ax.bar(df_mean.index, df_mean['NI'], color="#D0B9A4", width=0.6, linewidth=0.5, edgecolor='black', alpha=0.7)
    ax.errorbar(df_mean.index, df_mean['NI'], yerr=df_sd['NI'], fmt='none', ecolor='black', elinewidth=0.5, capsize=0)

    # 生データのドットプロット
    for AA in sample_name:
        available = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in df.index]
        x_positions = [AA] * len(available)

        y_values_NI = df.loc[available, "NI"].values#NIの生データ
        ax.scatter(x_positions, y_values_NI, color='black', s=10, edgecolor='black', linewidths=0.1, zorder=3, alpha=0.5)
  
    # NIの範囲設定
    ax.set_ylim(NI_YLIM)

    # 有意差のあるサンプルについて*マークをつける処理
    for i, (p_val, rect) in enumerate(zip(p_ni_adj, bars)):
        if p_val < 0.05:  # 有意差のある場合
            ax.text(rect.get_x() + rect.get_width() / 2., rect.get_y() - 3000, '*', ha='center', color='k', fontsize=10)
  
    # ラベルとタイトルの設定
    ax.set_title(f'NI: {TITLE}', fontsize=14)

    # プロットを表示
    #plt.show()
    plt.savefig(NIplot_file, dpi=350)

#=================================================================================================================================


Growth_df = tecan_data('Label1')
FI_df = tecan_data('Label2')
FI_df = subtract_phi(FI_df)

#Growth_df, FI_df = fourier_transform(Growth_df, GFP_df)
#------------------or------------------
create_lineplot(Growth_df, FI_df)
MGR, MFI, NI, MGR_pct, MFI_pct = calculation(Growth_df, FI_df)
plate_visualization(MGR, MFI, NI)
results = create_maindata(MGR, MFI, NI)
create_mainplot(results)
create_NIplot(results)