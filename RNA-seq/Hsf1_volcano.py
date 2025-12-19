import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats
import os

# Load the provided CSV file to check its contents
DEG_path = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig5/RNA-seq/kitamura_PolyI_reanalysis/Analysis/DEG_results.xlsx" #DEGデータを指定
gene_list = "/Users/muraseyukihiro/Desktop/PolyX_paper/new/fig5/RNA-seq/Gene_list_for_RNAseq.xlsx"
xls = pd.ExcelFile(DEG_path)
glist = pd.read_excel(gene_list)
path = os.path.dirname(DEG_path)

Hsf1 = ['AHA1', 'BTN2', 'CPR6', 'CUR1', 'FES1', 'HCH1', 'HSC82', 'HSP104', 
         'HSP42', 'HSP78', 'HSP82', 'MBF1', 'MDJ1', 'SIS1', 'SSA1', 'SSA2', 
         'STI1', 'YDJ1']

RPN4 = ['RPN4']

# Extract IDs for the specified genes
Hsf1_yname = glist[glist['Gene name'].isin(Hsf1)]['ID'].tolist()
RPN4_yname = glist[glist['Gene name'].isin(RPN4)]['ID'].tolist()


# volcano plotの作成--------------------------------------------------------------------------------
for sheet_name in xls.sheet_names:
    A, B = sheet_name.split("_vs_")
    df = pd.read_excel(xls, sheet_name=sheet_name)
    df['-log10(FDR)'] = -np.log10(df['FDR'])# -log10(FDR)を計算する

    df_Hsf1_points = df[df['Gene_id'].isin(Hsf1_yname)]
    #df_Msn24_points = df[df['Gene_id'].isin(Msn24_yname)]
    df_RPN4_points = df[df['Gene_id'].isin(RPN4_yname)]
    df_other_points = df[~df['Gene_id'].isin(Hsf1_yname + RPN4_yname)]

    # Volcano plotを描く 
    plt.figure(figsize=(8, 8))
    plt.scatter(df_other_points['logFC'], df_other_points['-log10(FDR)'], color='grey', alpha=0.3, edgecolors='none', label=None,s = 100)
    plt.scatter(df_Hsf1_points['logFC'], df_Hsf1_points['-log10(FDR)'], color="#FF0000", alpha=0.8, edgecolors='none', label=None,s = 150)
    plt.scatter(df_RPN4_points['logFC'], df_RPN4_points['-log10(FDR)'], color="#0000FF", alpha=0.8, edgecolors='none', label=None,s = 150)

    plt.axhline(-np.log10(0.05), color='gray', linestyle='--', label=None)# y=-log10(0.05) の線を追加

    #グラフデザイン
    ax = plt.gca()
    plt.xlim(-5, 5)
    plt.ylim(0, 60)

    # 原点クロス
    plt.axhline(y=0, color='black', linewidth=2)
    plt.axvline(x=0, color='black', linewidth=2)

    # 枠線調整
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.tick_params(width=2)

    # 目盛をゼロ軸に沿わせる
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # 文字サイズ
    #plt.xlabel("Log2FC", fontsize=15)
    #plt.ylabel("-log10(FDR)", fontsize=15)
    ax.tick_params(labelsize=15)

    plt.savefig(os.path.join(path, f'Stress_volcano_{A}_vs_{B}.png'), format='png', dpi=350)

#HSIxls.close()