import os
import pandas as pd
import time




def inx_count(ds):
    sum = ds.sum()
    sum_numerator = 0
    for count in ds:
        sum_numerator += count*(count -1)
    index = 1 - (sum_numerator / (sum *(sum -1)))
    return index

# path = r"D:\python\3.MLST_SC\第四次MLST\8.最终的MLST2\candidate_genes\candidate_genes.txt"
# genes = ["group"]
# with open(path) as file:
#     for line in file.readlines():
#         genes.append(line[0:-1])
# df_allels = pd.read_excel(r"D:\python\3.MLST_SC\第四次MLST\8.最终的MLST2\找最终组合\alles_strain.xlsx", index_col=0)
# df2 = df_allels[genes]
# print(df2)
# df2.to_excel("alles_strain_candidate_genes.xlsx")

df_group = pd.read_excel(r"D:\python\3.MLST_SC\第四次MLST\5.1找A克隆区分片段\6.candidate_genes\candidate_genes.xlsx")
df2 = df_group[["gene", "group"]]
groups_dict = dict()
for cha in list("01234567"):
    cha = int(cha)
    df3 = df2[df2.group == cha]
    groups_dict[cha] = list(df3.gene)
print(groups_dict)
groups_new = dict()
n =1
for gene1 in groups_dict[1]:
    for gene2 in groups_dict[2]:
        for gene3 in groups_dict[3]:
            for gene4 in groups_dict[4]:
                for gene5 in groups_dict[5]:
                    for gene6 in groups_dict[6]:
                        for gene7 in groups_dict[7]:
                            gene_group = [gene1, gene2, gene3, gene4, gene5, gene6, gene7, "group"]
                            groups_new[n] = gene_group
                            n += 1
print(groups_new)
print(len(groups_new))
print("各个组最小")
print("index最大", groups_new[112944])
#
#
df = pd.read_excel(r"D:\python\3.MLST_SC\第四次MLST\5.1找A克隆区分片段\6.candidate_genes\gene_allels.xlsx", index_col=0)


df_result = pd.DataFrame()
# n = 1
print(time.localtime(time.time()))
for k, group in groups_new.items():
    # n += 1
    # if n % 10000 == 0:
    #     print(n, time.localtime(time.time()))
    ds = pd.Series()
    df_group = df[group]
    # print(df_group)
    # 计算总的index
    df_all = df_group.drop(columns=["group"])
    index_all = inx_count(df_all.value_counts().values)
    ds["index_all"] = index_all
    # 计算各个group的index
    for cha in ["B", "C", "D", "E", "F", "L", "A", "A1", "A2", "A3", "A4", "A5", "A6", "Basal"]:
        df_g= df_group[df_group.group == cha].drop(columns=["group"])
        index = inx_count(df_g.value_counts().values)
        ds[cha] = index
    # 计算basal的index
    # df_basal = df_group[df_group.group == "Basal"].drop(columns=["group"])
    # index_basal = inx_count(df_basal.value_counts().values)
    # ds["basal"] = index_basal
    df_result[k] = ds
    # print(ds)

df2 = df_result.T
df2.to_excel("组合的index.xlsx")




