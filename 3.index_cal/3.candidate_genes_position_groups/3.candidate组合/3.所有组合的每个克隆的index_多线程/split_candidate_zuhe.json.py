# 读取所有candidate组合json文件，分为数个不同的csv文件

import os
import pandas as pd
import time
import ijson
from icecream import ic



# 读取所有的组合
path = r"../1.candidate_zuhe.json"
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

# for test
# gene_dict = {'1': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'alaS', 'czrB'],
#                 '10': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'secF', 'czrB'],
#                 '11': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'tgt', 'czrB'],
#                 '12': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'thiI', 'czrB'],
#                 '13': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'tig', 'czrB'],
#                 '14': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'alaS', 'czrB'],
#                 '15': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'dnaI', 'czrB'],
#                 '16': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'gatA', 'czrB'],
#                 '17': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'groEL', 'czrB'],
#                 '18': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'metK', 'czrB'],
#                 '19': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'murC', 'czrB'],
#                 '2': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'dnaI', 'czrB'],
#                 '20': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'obgE', 'czrB'],
#                 '21': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'rarA', 'czrB'],
#                 '22': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'scrB', 'czrB'],
#                 '23': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'secF', 'czrB'],
#                 '24': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'tgt', 'czrB'],
#                 '25': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'thiI', 'czrB'],
#                 '26': ['acsA2', 'guaB', 'addB', 'carB', 'aroB', 'tig', 'czrB'],
#                 '27': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'alaS', 'czrB'],
#                 '28': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'dnaI', 'czrB'],
#                 '29': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'gatA', 'czrB'],
#                 '3': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'gatA', 'czrB'],
#                 '30': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'groEL', 'czrB'],
#                 '31': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'metK', 'czrB'],
#                 '32': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'murC', 'czrB'],
#                 '33': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'obgE', 'czrB'],
#                 '34': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'rarA', 'czrB'],
#                 '35': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'scrB', 'czrB'],
#                 '36': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'secF', 'czrB'],
#                 '37': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'tgt', 'czrB'],
#                 '38': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'thiI', 'czrB'],
#                 '39': ['acsA2', 'guaB', 'addB', 'carB', 'birA', 'tig', 'czrB'],
#                 '4': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'groEL', 'czrB'],
#                 '40': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'alaS', 'czrB'],
#                 '41': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'dnaI', 'czrB'],
#                 '42': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'gatA', 'czrB'],
#                 '43': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'groEL', 'czrB'],
#                 '44': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'metK', 'czrB'],
#                 '45': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'murC', 'czrB'],
#                 '46': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'obgE', 'czrB'],
#                 '47': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'rarA', 'czrB'],
#                 '48': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'scrB', 'czrB'],
#                 '49': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'secF', 'czrB'],
#                 '5': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'metK', 'czrB'],
#                 '50': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'tgt', 'czrB'],
#                 '51': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'thiI', 'czrB'],
#                 '52': ['acsA2', 'guaB', 'addB', 'carB', 'dinG', 'tig', 'czrB'],
#                 '53': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'alaS', 'czrB'],
#                 '54': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'dnaI', 'czrB'],
#                 '55': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'gatA', 'czrB'],
#                 '56': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'groEL', 'czrB'],
#                 '57': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'metK', 'czrB'],
#                 '58': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'murC', 'czrB'],
#                 '59': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'obgE', 'czrB'],
#                 '6': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'murC', 'czrB'],
#                 '60': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'rarA', 'czrB'],
#                 '61': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'scrB', 'czrB'],
#                 '62': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'secF', 'czrB'],
#                 '63': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'tgt', 'czrB'],
#                 '64': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'thiI', 'czrB'],
#                 '65': ['acsA2', 'guaB', 'addB', 'carB', 'lysC', 'tig', 'czrB'],
#                 '66': ['acsA2', 'guaB', 'addB', 'carB', 'odhB', 'alaS', 'atpB_2'],
#                 '67': ['acsA2', 'guaB', 'addB', 'carB', 'odhB', 'alaS', 'cbiO_1'],
#                 '68': ['acsA2', 'guaB', 'addB', 'carB', 'odhB', 'alaS', 'ctrA'],
#                 '7': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'obgE', 'czrB'],
#                 '8': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'rarA', 'czrB'],
#                 '9': ['acsA2', 'guaB', 'addB', 'carB', 'aroA', 'scrB', 'czrB']}

df = pd.DataFrame(gene_dict).T
line_num = df.shape[0]
output_dir = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\3.candidate组合\3.所有组合的每个克隆的index_多线程\split_json"
split_num = 12
group_length = line_num // split_num
if line_num % split_num != 0:
    group_length += 1
ic(group_length)
for a in range(split_num):
    file_name = f"candidate_zuhe_{a}.csv"
    file_path = os.path.join(output_dir, file_name)
    start = a * group_length
    end = (a+1) * group_length
    df_2 = df.iloc[start:end]
    df_2.to_csv(file_path)


