import pandas as pd
import numpy as np


strains = ["BY4741", "cdc24", "mac1", "mmr1", "pre7", "psk1", "rpl18b", "rpl19a"]
plasmids = ["phi", "delta", "PolyE", "mox"]
output_file = '/Users/muraseyukihiro/Desktop/PolyX/Microscope/elongation_PolyE/all_data/count_result.csv'  # 保存するファイル名

df = pd.read_csv("/Users/muraseyukihiro/Desktop/PolyX/Microscope/elongation_PolyE/all_data/result.csv")

result_df = pd.DataFrame(columns=plasmids, index=strains)

#変異株ごとに数える
for strain in strains:
    for plasmid in plasmids:
        count = len(df[(df['strain'] == strain) & (df['plasmid'] == plasmid)])
        count_long = len(df[(df['strain'] == strain) & (df['plasmid'] == plasmid) & (df['Major/Minor'] >= 1.5)])

        ratio = count_long / count
        result_df.loc[strain, plasmid] = ratio

result_df.to_csv(output_file)