import os
import pandas as pd
import matplotlib.pyplot as plt

path_dir = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\1.best_fragment_index_json\2.useful_csv"
result = {"gene": [], "fragments_num": []}

for file in os.listdir(path_dir):
    gene = file.split(".")[0]
    path_csv = os.path.join(path_dir, file)
    # print(path_csv)
    df = pd.read_csv(path_csv)
    result["gene"].append(gene)
    result["fragments_num"].append(df.shape[0])
print(result)
df2 = pd.DataFrame(result)
# df2.to_excel("useful_csv_stastic.xlsx")
print(df2.describe())

a = df2.fragments_num.value_counts()
print(a)
print(a.describe())