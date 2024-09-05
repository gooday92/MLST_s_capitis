import multiprocessing
import pandas as pd
from icecream import ic
import json
import os


def inx_count(ds):
    sum = ds.sum()
    sum_numerator = 0
    for count in ds:
        sum_numerator += count*(count -1)
    index = 1 - (sum_numerator / (sum *(sum -1)))
    return index

def groups_zuhe(df):
    groups = []
    for index, row in df.iterrows():
        group = list(row)
        group.append("group")
        groups.append(group)
    return groups


def index_all(df_alleles, df_groups, out_file):
    # 输入 alleles的df，组合的df,
    # 输出 不同组合的不同克隆的index,保存为csv
    # ic(df_alleles, df_groups)
    # ic(df_alleles, df_groups)
    # 遍历group组合
    df_result = pd.DataFrame()
    groups = groups_zuhe(df_groups)
    for group in groups:
        # quality = 1
        k = "".join(group)
        ds = dict()
        df_sub = df_alleles[group]
        # 计算总的index
        df_all = df_sub.drop(columns=["group"])
        index_all = inx_count(df_all.value_counts().values)
        ds["index_all"] = index_all
        # 计算各个clone的index
        for cha in ["B", "C", "D", "E", "F", "L", "A", "Basal"]:
            df_g = df_sub[df_sub.group == cha].drop(columns=["group"])
            index = inx_count(df_g.value_counts().values)
            # if index > 0.1:
            #     quality = 0
            ds[cha] = index
        # if quality == 1:
        df_result[k] = pd.Series(ds)
    df_result = df_result.T
    df_result.to_csv(out_file)



if __name__ == "__main__":
    json_dir = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\3.candidate组合\3.所有组合的每个克隆的index_多线程\split_json"
    out_dir = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\3.candidate组合\3.所有组合的每个克隆的index_多线程\index_out_csv"
    split_num = 12
    df_alleles = pd.read_excel(r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\gene_allels.xlsx", index_col=0)
    pool = multiprocessing.Pool(split_num)
    for a in range(split_num):
        file_name = f"candidate_zuhe_{a}.csv"
        file_path = os.path.join(json_dir, file_name)
        out_path = os.path.join(out_dir, file_name)
        df = pd.read_csv(file_path, index_col=0)
        pool.apply_async(func=index_all, args=(df_alleles, df, out_path))
    pool.close()
    pool.join()
