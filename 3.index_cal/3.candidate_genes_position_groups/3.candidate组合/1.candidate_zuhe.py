import json
import os
import pandas as pd
import time
from icecream import ic

df2 = pd.read_excel(r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\2.position\position_groups_result.xlsx")
# ic(df2)
groups_dict = dict()
for cha in list("01234567"):
    cha = int(cha)
    df3 = df2[df2.position == cha]
    groups_dict[cha] = list(df3.name)
print(groups_dict)

df2.set_index("name", drop=False, inplace=True)
ic(df2)
groups_new = dict()
n = 1
for gene1 in groups_dict[1]:
    for gene2 in groups_dict[2]:
        for gene3 in groups_dict[3]:
            for gene4 in groups_dict[4]:
                for gene5 in groups_dict[5]:
                    for gene6 in groups_dict[6]:
                        for gene7 in groups_dict[7]:
                            gene_group = [gene1, gene2, gene3, gene4, gene5, gene6, gene7]
                            df = df2.loc[gene_group, :]
                            a = ""
                            for l in df.clones.unique():
                                a += l
                            b = set(list(a))
                            if len(b) == 7:
                                groups_new[n] = gene_group
                                n += 1
# print(groups_new)
print(len(groups_new))
a = json.dumps(groups_new,indent=2)
f2 = open("candidate_zuhe.json", "w")
f2.write(a)
f2.close()
