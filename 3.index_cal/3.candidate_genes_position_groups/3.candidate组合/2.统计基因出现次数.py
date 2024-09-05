import json
import pandas as pd
from icecream import ic


path = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\3.candidate组合\1.candidate_zuhe.json"
file = open(path)
f = json.load(file)
file.close()
df = pd.DataFrame(f)
df.to_csv("1.candidate_zuhe.csv")
ic(f)
ic(df)


