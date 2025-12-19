import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
plt.rcParams["font.family"] = "Arial"
import openpyxl
from ipywidgets import interact, FloatSlider

max_df = pd.read_excel("/Users/muraseyukihiro/Desktop/PolyX/Dry/makino/various_polyX.xlsx", sheet_name="max_polyX", index_col=0)
num_df = pd.read_excel("/Users/muraseyukihiro/Desktop/PolyX/Dry/makino/various_polyX.xlsx", sheet_name="num_poly10X per gene", index_col=0)
num_df = num_df*1000

def show_pal2(start, rot):
    sns.palplot(sns.cubehelix_palette(24, start=start, rot=rot))
interact(show_pal2, start=FloatSlider(max=1), rot=FloatSlider(0.4, min=-1, max=1))

sns.choose_colorbrewer_palette('diverging')

A = sns.color_palette("RdBu_r", 20, desat=0.7)
C = sns.color_palette("PiYG_r", 3, desat=0)
D = sns.color_palette("viridis", 20, desat=0.7)

#max_polyX解析

plt.figure(dpi = 300, figsize=(21,15))

sns.set(font_scale=1.9)

ax = sns.heatmap(max_df,
                vmin=8, vmax=40,
                annot=True, fmt='d',
                cmap=A)

plt.yticks(rotation=360)
plt.xticks(rotation=30)

plt.savefig("/Users/muraseyukihiro/Desktop/PolyX/Dry/makino/max_polyX.pdf")

# num_poly10X解析

#マスクを作成する
num_df_mask = (num_df == 0)

#ここからfigure
plt.figure(dpi = 300, figsize=(21,15))

sns.set(font_scale=1.9)

ax = sns.heatmap(num_df,
                vmin=0, vmax=4,
                annot=True,
                cmap=C, 
                cbar=False)

#maskを適用（Falseの場所に色を入れる）
sns.heatmap(num_df,
            vmin=0, vmax=4,
            annot=True,
            cmap=A, 
            mask=num_df_mask )

plt.yticks(rotation=360)
plt.xticks(rotation=30)

plt.savefig("/Users/muraseyukihiro/Desktop/PolyX/Dry/makino/num_poly10X.pdf")
