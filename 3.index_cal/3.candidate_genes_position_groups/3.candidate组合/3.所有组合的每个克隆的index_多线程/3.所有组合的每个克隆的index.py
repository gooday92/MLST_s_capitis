import os
import pandas as pd
import time
import ijson
from icecream import ic

def inx_count(ds):
    sum = ds.sum()
    sum_numerator = 0
    for count in ds:
        sum_numerator += count * (count - 1)
    index = 1 - (sum_numerator / (sum * (sum - 1)))
    return index

# 读取所有的组合
path = r"candidate_zuhe.json"
file = open(path)
data = ijson.parse(file)
gene_dict = dict()
id = ""
for num, type, gene in data:
    if type == 'string':
        gene_dict[id].append(gene)
    elif type == "start_array":
        id = num
        gene_dict[id] = []

df = pd.read_excel(r"gene_allels.xlsx", index_col=0)
df_result = pd.DataFrame()
