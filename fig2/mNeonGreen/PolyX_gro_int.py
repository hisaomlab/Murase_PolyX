import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smm
plt.rcParams["font.family"] = "Arial"

#外れデータはサンプルシートでEmptyとして指定する

#-----読み込むファイルの指定-----
TECAN_file_name = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/mNeonGreen/20241217_mNeonGreen_PolyX.xlsx" #TECANのraw_dataを指定
sample_sheet_name = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/mNeonGreen/Sample_Sheet_Poly10X.xlsx" #サンプルシートを指定
#sample_name = ['Δ','D', 'Q', 'N', 'S', 'E', 'R', 'H', 'P', 'K', 'A', 'V', 'T','G', 'L', 'C','F','I','Y','M','W','Φ'] #サンプル名前を指定
sample_name = ['Δ','A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N','P', 'Q', 'R','S','T','V','W','Y'] #サンプル名前を指定
path = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig2/mNeonGreen"#保存先のパス
plate_visualization_file = f'{path}/platedata.xlsx' #最大蛍光値・MGRを保存するファイル
lineplot_file = f'{path}/lineplot.pdf' #lineplotを保存するファイル
mainplot_file = f'{path}/mainplot.pdf' #mainplotを保存するファイル
NIplot_file = f'{path}/NIplot.pdf'#NIplotを保存するファイル
miandata_file = f'{path}/maindata.xlsx'#maindataを保存するファイル
#-----------------------------

TECAN_df = pd.read_excel(TECAN_file_name, sheet_name=0, header=None)

#サンプルシートを読み込み、wellとsample_numberを対応させる
sample_sheet_df = pd.read_excel(sample_sheet_name, sheet_name=0)
well_to_sample = dict(zip(sample_sheet_df['Well'], sample_sheet_df['Sample_number']))

#TECANのデータから各項目の範囲を読み込み、Sample_number列と結合する--------------------------------------------
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
#------------------------------------------------------------------------------------------------------

def fourier_transform(Growth_df, GFP_df):
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
    GFP_ave_df =  GFP_df.rolling(period,axis=1,center=True,min_periods=1).mean()

    create_lineplot(Growth_ave_df, GFP_ave_df)
    calculation(Growth_ave_df, GFP_ave_df)
    
#増殖、蛍光曲線をグラフで表す============================================================
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

def create_lineplot(growth_df, gfp_df):
    fin_time = 80  # [hour]
    fig = plt.figure(dpi=100, figsize=(20, 10))
    time = np.arange(len(growth_df.T)) / 6  # interval time = 10min の場合
    #time = np.arange(len(growth_df.T)) / 2  # interval time = 30min の場合

    for n, AA in zip(range(len(sample_name)), sample_name):
        ax = fig.add_subplot(4, 6, n + 1)
        available_samples = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in growth_df.index]#_1,_2,_3,_4の中であるものを取得する

        # 増殖曲線グラフ
        growth_data = growth_df.loc[available_samples]
        plot_filled_line(ax, AA, growth_data, time, color="#23375B", label="_nolegend_", ylim=(0, 1.5), yticks=[0.5, 1.0, 1.5])

        # GFP蛍光強度のグラフ
        ax2 = ax.twinx()
        gfp_data = gfp_df.loc[available_samples]# / 10  # 蛍光値を10分の1にスケール
        plot_filled_line(ax2, AA, gfp_data, time, color="#58C9A7", label="_nolegend_", ylim=(0, 15000), yticks=[5000, 10000, 15000])

        # x軸の表示範囲
        ax.set_xlim(0, fin_time)
        ax.set_xticks([0, fin_time / 2, fin_time])

        # グラフのスタイル設定
        for axis in [ax, ax2]:
            axis.spines["top"].set_visible(False)
            axis.spines["left"].set_linewidth(2)
            axis.spines["right"].set_linewidth(2)
            axis.spines["bottom"].set_linewidth(2)
            axis.tick_params(width=1, labelsize=12)

    plt.subplots_adjust(wspace=0.45, hspace=0.45)
    #plt.show()
    plt.savefig(lineplot_file, dpi=350)

#==========================================================================================

def calculation(Growth_df, GFP_df):
    #~~~~パラメータ~~~~
    Period = 10 #測定間隔（min）
    Range = 25 #GRを算出するレンジ
    Step_size = 10 #移動測定数（大きくすると処理が早くなる）
    Threshold = 0.9 #R2の閾値を指定
    #MGR算出------------------------------------------------------------------------------------------------------
    #データをログ変換
    OD_log2 = np.log2(Growth_df).dropna()
    Max_cycle = len(OD_log2.columns)
    Well_Num = len(OD_log2)
    x = np.arange(Range)

    MGR = OD_log2[1]

    #WellごとにTimeウインドウを移動しながらGR（直線回帰の傾き=SlopeとR2を取得）、R2のThreshold以上のSlopeを返す（Threshold以下だと0.01を返す）。
    #最終的なGRは測定時間間隔（Period）で割り(min-1)、便宜的に1000をかけている。

    for index in MGR.index:
        Cycle = 70 #GR計算開始の位置(計算を始めたいところのサイクル数を入れる)
        n = 1
        Slope  = np.array([0.1])
        R2 = np.array([1.0])
        OD_list = OD_log2.loc[index].values

        while Cycle < Max_cycle-Range:
            OD_list2 = OD_list[Cycle:Cycle+Range]
            Slope = np.insert(Slope,n,(np.polyfit(x,OD_list2,1)[0] * 1000 / Period)) #1次関数で傾き（Slope）を計算。1000倍してPeriod(min)で割り、Slopeに順次代入
            R2 = np.insert(R2,n,(np.corrcoef(x,OD_list2)[0][1]**2)) #R2を計算。R2に順次代入
            Cycle += Step_size
            n += 1

        GR = np.array([x for x in Slope[R2 > Threshold]]) #R2がThresholdより高いGRのarrayを作成
        MGR[index] = GR.max()

    GFP = GFP_df.max(axis=1)#最大GFP蛍光強度

    NI = MGR*GFP#中性指数

    plate_visualization(MGR, GFP, NI)
    create_maindata(MGR, GFP, NI)

#-----------------------------------------------------------------------------------------------------

def plate_visualization(MGR, GFP, NI):
    #MGRの整形
    MGR = np.array(MGR.values)
    MGR = MGR.reshape([8,12])

    #GFPの整形
    GFP = np.array(GFP.values)
    GFP = GFP.reshape([8,12])

    #NIの整形
    NI = np.array(NI.values)
    NI = NI.reshape([8,12])

    column = [i for i in range(1, 13)]
    row = ["A", "B", "C", "D", "E", "F", "G", "H"]

    MGR = pd.DataFrame(MGR,columns=column, index=row)
    GFP = pd.DataFrame(GFP,columns=column, index=row)
    NI = pd.DataFrame(NI,columns=column, index=row)

    with pd.ExcelWriter(plate_visualization_file) as writer:
        MGR.to_excel(writer,sheet_name="MGR")
        GFP.to_excel(writer,sheet_name="GFP")
        NI.to_excel(writer,sheet_name="NI(MGR_GFP)")

def create_maindata(MGR, GFP, NI):
    df = pd.DataFrame({"MGR": MGR, "GFP": GFP, "NI": NI})
    df_mean = pd.DataFrame(columns=["MGR", "GFP", "NI"])
    df_sd = pd.DataFrame(columns=["MGR", "GFP", "NI"])
    df_p = pd.DataFrame(columns=["MGR", "GFP", "NI"])

    # p値を出す際の基準(mox)のサンプルを取得
    base_samples = [s for s in ["Δ_1", "Δ_2", "Δ_3", "Δ_4"] if s in df.index]

    for AA in sample_name:
        available_samples = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in df.index]#_1,_2,_3,_4の中であるものを取得する
        MGR_mean = MGR.loc[available_samples].mean()
        GFP_mean = GFP.loc[available_samples].mean()
        NI_mean = NI.loc[available_samples].mean()
        df_mean.loc[AA] = {"MGR": MGR_mean, "GFP": GFP_mean, "NI": NI_mean}

        MGR_sd = MGR.loc[available_samples].std()
        GFP_sd = GFP.loc[available_samples].std()
        NI_sd = NI.loc[available_samples].std()
        df_sd.loc[AA] = {"MGR": MGR_sd, "GFP": GFP_sd, "NI": NI_sd}

        t_stat, p_value_mgr = ttest_ind(MGR.loc[available_samples], MGR.loc[base_samples], equal_var=True)# moxとの比較MGR (t検定)
        t_stat, p_value_gfp = ttest_ind(GFP.loc[available_samples], GFP.loc[base_samples], equal_var=True)# moxとの比較GFP (t検定)
        t_stat, p_value_ni = ttest_ind(NI.loc[available_samples], NI.loc[base_samples], equal_var=True)# moxとの比較NI (t検定)
        df_p.loc[AA] = {"MGR" : p_value_mgr, "GFP" : p_value_gfp, "NI" : p_value_ni}

    df_mean["GFP_Φ"] = (df_mean["GFP"] - df_mean.loc["Φ", "GFP"]).apply(lambda x: max(x, 0))#Φ(自家蛍光)分を排除する
    
    # ボンフェロー二補正
    p_values_mgr = df_p["MGR"].values
    p_values_gfp = df_p["GFP"].values
    p_values_ni = df_p["NI"].values
    _, p_adj_mgr, _, _ = smm.multipletests(p_values_mgr, alpha=0.05, method='bonferroni')
    _, p_adj_gfp, _, _ = smm.multipletests(p_values_gfp, alpha=0.05, method='bonferroni')
    _, p_adj_ni, _, _ = smm.multipletests(p_values_ni, alpha=0.05, method='bonferroni')
    # ボンフェローニ補正結果をデータフレームにまとめる
    df_p["MGR_bonferroni"] = p_adj_mgr
    df_p["GFP_bonferroni"] = p_adj_gfp
    df_p["NI_bonferroni"] = p_adj_ni

    # Excelにデータフレームを書き込む
    with pd.ExcelWriter(miandata_file) as writer:
        df.to_excel(writer, sheet_name="Data")
        df_mean.to_excel(writer, sheet_name="Mean")
        df_sd.to_excel(writer, sheet_name="SD")
        df_p.to_excel(writer, sheet_name="P-Values")

    create_mainplot(MGR, GFP, df, df_mean, df_sd, p_adj_mgr, p_adj_gfp)
    create_NIplot(NI, df, df_mean, df_sd, p_adj_ni)

#MGR,MFI,NIを棒グラフで表す=======================================================================
def create_mainplot(MGR, GFP, df, df_mean, df_sd, p_adj_mgr, p_adj_gfp):
    fig, ax = plt.subplots(figsize=(10, 3))# プロットの準備

   # MGR 平均値と標準偏差の棒グラフ
    bars = ax.bar(df_mean.index, df_mean['MGR'], color=np.where(df_mean.index.isin(["Φ"]), "gray", "#23375B"), label='Max growth rate', width=0.6, linewidth=0.5, edgecolor='black', zorder=1, alpha=0.7)
    ax.errorbar(df_mean.index, df_mean['MGR'], yerr=df_sd['MGR'], fmt='none', ecolor='black', elinewidth=0.5, capsize=0, zorder=2)

    # GFP 平均値と標準偏差のプロット
    ax2 = ax.twinx()
    ax2.scatter(df_mean.index, df_mean['GFP_Φ'], color='#58C9A7', marker='_', label='Expression level', s=300, edgecolor='black', linewidth=3, zorder=2, alpha=0.9)
    ax2.errorbar(df_mean.index, df_mean['GFP_Φ'], yerr=df_sd['GFP'], fmt='none', ecolor='black', elinewidth=0.5, capsize=0, zorder=1)

    # 生データのドットプロット
    for AA in sample_name:
        available_samples = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in df.index]
        x_positions = [AA] * len(available_samples)

        y_values_MGR = MGR.loc[available_samples].values#MGRの生データ
        ax.scatter(x_positions, y_values_MGR, color='black', s=10, linewidths=0.1, zorder=3, alpha=0.5)

        y_values_GFP = GFP.loc[available_samples].values#GFPの生データ
        y_values_GFP = np.maximum(y_values_GFP - df_mean.loc["Φ", "GFP"], 0)
        ax2.scatter(x_positions, y_values_GFP, color='#58C9A7', s=10, linewidths=0.1, zorder=3, alpha=0.7)
  
    # MGRの範囲設定
    ax.set_ylim([0, 5])  # MGRの範囲で設定
    ax2.set_ylim([0, 15000])  # GFPの範囲で設定

    # 有意差のあるサンプルについて*マークをつける処理(MGR)
    for i, (p_val, rect) in enumerate(zip(p_adj_mgr, bars)):
        if p_val < 0.05:  # 有意差のある場合
            ax.text(rect.get_x() + rect.get_width() / 2., rect.get_y() - 0.75, '*', ha='center', color='k', fontsize=10)
    
    # 有意差のあるサンプルについて*マークをつける処理(GFP)
    for i, (p_val, rect) in enumerate(zip(p_adj_gfp, bars)):
        if p_val < 0.05:  # 有意差のある場合
            ax.text(rect.get_x() + rect.get_width() / 2., rect.get_y() - 0.8, '*', ha='center', color='green', fontsize=10)
  
    # ラベルとタイトルの設定
    #ax.set_ylabel('Max growth rate', fontsize=12)
    #ax2.set_ylabel('Expression level', fontsize=12, rotation=270)
    ax.set_title('MGR and MFI: mNeonGreen-PolyX (SC-LU)', fontsize=14)

    # プロットを表示
    #plt.show()
    plt.savefig(mainplot_file, dpi=350)

def create_NIplot(NI, df, df_mean, df_sd, p_adj_ni):
    fig, ax = plt.subplots(figsize=(10, 3))# プロットの準備

   # NI 平均値と標準偏差の棒グラフ
    bars = ax.bar(df_mean.index, df_mean['NI'], color="#D0B9A4", width=0.6, linewidth=0.5, edgecolor='black', alpha=0.7)
    ax.errorbar(df_mean.index, df_mean['NI'], yerr=df_sd['NI'], fmt='none', ecolor='black', elinewidth=0.5, capsize=0)

    # 生データのドットプロット
    for AA in sample_name:
        available_samples = [s for s in [AA + "_1", AA + "_2", AA + "_3", AA + "_4"] if s in df.index]
        x_positions = [AA] * len(available_samples)

        y_values_NI = NI.loc[available_samples].values#NIの生データ
        ax.scatter(x_positions, y_values_NI, color='black', s=10, edgecolor='black', linewidths=0.1, zorder=3, alpha=0.5)
  
    # MGRの範囲設定
    ax.set_ylim([0, 50000])  # MGRの範囲で設定

    # 有意差のあるサンプルについて*マークをつける処理
    #for i, (p_val, rect) in enumerate(zip(p_adj_ni, bars)):
    #    if p_val < 0.05:  # 有意差のある場合
    #        ax.text(rect.get_x() + rect.get_width() / 2., rect.get_y() - 110, '*', ha='center', color='k', fontsize=16)
  
    # ラベルとタイトルの設定
    #ax.set_ylabel('Neutral Index (MGR*GFP)', fontsize=12)
    ax.set_title('NI: mNeonGreen-PolyX (SC-LU)', fontsize=14)

    # プロットを表示
    #plt.show()
    plt.savefig(NIplot_file, dpi=350)

#=================================================================================================================================


Growth_df = tecan_data('Label1')
GFP_df = tecan_data('Label2')

#fourier_transform(Growth_df, GFP_df)
#------------------or------------------
create_lineplot(Growth_df, GFP_df)
calculation(Growth_df, GFP_df)