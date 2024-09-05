import os
import pandas as pd
from icecream import ic

gene_list = ['phoA', 'mntC', 'clpP', 'carB', 'rluB', 'hisS', 'atpB_2']

df = pd.read_excel("gene_allels.xlsx", index_col=0)
df = df[gene_list]

strain_st = dict()
for starin, row in df.iterrows():
    # print(row)
    st = ""
    for i in row.values:
        st += str(i)
    strain_st[starin] = st
ds = pd.Series(strain_st)
l = list(ds.value_counts().index)
st_types = dict()
for i, st in enumerate(l):
    st_types[st] = i+1

strain_types = dict()
for strain, st in strain_st.items():
    strain_types[strain] = st_types[st]

df["ST"] = pd.Series(strain_types)
df.to_excel("strain_ST.xlsx")
ic(df)


