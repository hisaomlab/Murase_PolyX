import os
import gzip
import re
import pandas as pd
from Bio.PDB import PDBParser
from collections import defaultdict

# === 設定 ===
target_folder = "/Users/muraseyukihiro/Desktop/PolyX/Dry/Where_is_PolyX_in_protein/UP000002311_559292_YEAST_v4"
output = "/Users/muraseyukihiro/Desktop/PolyX/Dry/Where_is_PolyX_in_protein/Where_is_polyX_results.csv"
min_repeat = 1
max_repeat = 50

residue_3to1 = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
    'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',
    'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
    'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'
}

all_results = []

for filename in os.listdir(target_folder):
    if not filename.endswith(".pdb.gz"):
        continue

    filepath = os.path.join(target_folder, filename)
    parser = PDBParser(QUIET=True)

    try:
        with gzip.open(filepath, "rt") as f:
            structure = parser.get_structure("model", f)
    except Exception as e:
        print(f"Failed to read {filename}: {e}")
        continue

    sequence = ""
    residue_ids = []
    residue_to_bfactors = defaultdict(list)

    try:
        model = list(structure)[0]
        chain = list(model)[0]
    except:
        print(f"No valid model/chain in {filename}")
        continue

    for residue in chain:
        resname = residue.get_resname()
        if resname in residue_3to1:
            aa = residue_3to1[resname]
            sequence += aa
            res_id = residue.get_id()[1]
            residue_ids.append(res_id)
            bfactors = [atom.get_bfactor() for atom in residue if atom.get_name() == "CA"]
            residue_to_bfactors[res_id] = bfactors

    used_positions = set()  # ← すでに使ったインデックスを記録

    # 長い繰り返しから優先的に処理（重複排除のため）
    for repeat_len in range(max_repeat, min_repeat - 1, -1):
        pattern = '|'.join([f"{aa}{{{repeat_len},}}" for aa in "ACDEFGHIKLMNPQRSTVWY"])
        matches = re.finditer(pattern, sequence)

        for match in matches:
            start = match.start()
            end = match.end()

            # この範囲が既に使われていたらスキップ
            if any(i in used_positions for i in range(start, end)):
                continue

            poly_seq = match.group()
            poly_type = poly_seq[0]

            start_res = residue_ids[start]
            end_res = residue_ids[end - 1]

            scores = []
            for i in range(start, end):
                used_positions.add(i)  # 今後この範囲を除外
                scores.extend(residue_to_bfactors[residue_ids[i]])

            if scores:
                avg_score = sum(scores) / len(scores)
                max_score = max(scores)
                min_score = min(scores)
            else:
                avg_score = max_score = min_score = None

            all_results.append({
                "File": filename,
                "PolyX Type": poly_type,
                "Repeat Length": repeat_len,
                "Start Residue": start_res,
                "End Residue": end_res,
                "Length": len(poly_seq),
                "Avg pLDDT": avg_score,
                "Max pLDDT": max_score,
                "Min pLDDT": min_score
            })
            print(f"finish {filename}")

# Excel出力
df = pd.DataFrame(all_results)
df.sort_values(by=["File", "Repeat Length", "Start Residue"], inplace=True)
df.to_csv(output, index=False)
print(f"結果を {output} に保存しました。")
