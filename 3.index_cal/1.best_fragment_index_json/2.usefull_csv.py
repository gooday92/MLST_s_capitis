import os
import pandas as pd
import shutil


dir_from = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\1.best_fragment_index_json\1.slid_json"
dir_to = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\1.best_fragment_index_json\2.useful_json"


useful_list = []
for f in os.listdir(dir_from):
    gene = f.split(".")[0]
    file_path = os.path.join(dir_from, f)
    df = pd.read_csv(file_path)
    if len(df.columns) > 1:
        useful_list.append(f)
print(len(useful_list))
for f in useful_list:
    source = os.path.join(dir_from, f)
    target = os.path.join(dir_to, f)
    shutil.copy(source, target)




