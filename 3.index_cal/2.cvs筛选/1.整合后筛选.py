import os
import pandas as pd
from icecream import ic

dir_csv = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\1.best_fragment_index_json\2.useful_csv"

df2 = pd.DataFrame()
for f in os.listdir(dir_csv):
    gene_name = f.split(".")[0]
    file_path = os.path.join(dir_csv, f)
    df = pd.read_csv(file_path)
    df["gene"] = gene_name
    df.rename(columns={"Unnamed: 0": "start"}, inplace=True)
    df2 = pd.concat([df2, df], axis=0)
ic(df2)
df2 = df2[(df2.A < 0.1) & (df2.B < 0.1) & (df2.C < 0.1) & (df2.D < 0.1) & (df2.E < 0.1) & (df2.F < 0.1) & (df2.L < 0.1) & (df2.All > 0.3)]
# 143个符合标准的基因
ic(df2)
ic(df2.gene.unique())
ic(df2.gene.value_counts())
df2.to_excel("1_primary_fileter.xlsx", index=False)

