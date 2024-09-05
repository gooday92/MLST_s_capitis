import os
import pandas as pd
from icecream import ic

# function
def allels(a):
    return a.split("(")[0]



file = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\2.cvs筛选\1_primary_fileter.xlsx"
df = pd.read_excel(file)
num_allels_dict = dict()
for index, row in df.iterrows():
    l = ["A_num", "B_num", "C_num", "D_num", "E_num", "F_num", "L_num"]
    num_allels = len(row[l].map(allels).unique())
    num_allels_dict[index] = num_allels
    # map function
ds = pd.Series(num_allels_dict)
# ic(ds)
df["alleles_num"] = ds
df = df[df.alleles_num >= 4]
# 剩下136个基因，1534个片段
# ic(df)
# ic(df.gene.unique())
# ic(len(df.gene.unique()))
# ic(df.gene.value_counts())
# df.to_excel("2.allels_num.xlsx", index=False)

# # 每个基因内选择等位基因数目最高的片段
df_n = pd.DataFrame()
for gene in df.gene.unique():
    df2 = df[df.gene == gene]
    df2["sum"] = df.A + df.B + df.C + df.D + df.E + df.F + df.L
    # 每个基因取各个克隆的综合index最大
    df2.sort_values(by=["alleles_num", "sum"], ascending=[False, False])
    ds = df2.iloc[0, :]
    gene, start = ds["gene"], ds["start"]
    index = f"{gene}_{start}"
    df_n[index] = ds

df_n = df_n.T
df_n = df_n[df_n["sum"] < 0.15]
ic(df_n)
ic(df_n.gene.unique())
ic(len(df_n.gene.unique()))
ic(df_n.gene.value_counts())
# 剩下30个基因，30个片段
df_n.to_excel("2.3_candidate_genes.xlsx", index=False)
# #
clones_dict = dict()
for index, row in df_n.iterrows():
    l = ["A_num", "B_num", "C_num", "D_num", "E_num", "F_num", "L_num"]
    clones = list(row[l].map(allels).drop_duplicates(keep=False).index.map(lambda a: a[0]))
    c = ""
    for a in clones:
        c += a
    clones_dict[index] = c
ds = pd.Series(clones_dict)
df_n["clones"] = ds
ic(df_n)
df_n.to_excel("2.4_candidate_genes_clones.xlsx", index=False)
