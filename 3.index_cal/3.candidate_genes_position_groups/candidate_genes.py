from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import json
from icecream import ic


def align_df(path):
    alignment = AlignIO.read(path, "fasta")
    genes_list = []
    for record in alignment:
        genes_list.append(record.id)
    df = pd.DataFrame(alignment, index=genes_list)
    return df


def allel(df):
    type_dict = dict()
    for strain, row in df.iterrows():
        type_seq = "".join(row)
        strain = strain.split(";")[0]
        type_dict[strain] = type_seq
    ds = pd.Series(type_dict)
    df_groups = pd.DataFrame()
    df_groups["type"] = ds
    # ic(df_groups)
    # 编号不同的等位基因
    l = list(df_groups.type.value_counts().index)
    allels_type = dict()
    for i, seq in enumerate(l):
        allels_type[seq] = i + 1
    allels_genes = dict()
    for gene, row in df_groups.iterrows():
        allels_genes[gene] = allels_type[row.type]
    ds = pd.Series(allels_genes)
    return ds


def consensus_seq(df):
    seq = []
    for label, col in df.iteritems():
        # 取每一列占比最大的那个碱基
        cha = col.value_counts().index[0]
        seq.append(cha)
    conses_seq = "".join(seq)
    return conses_seq


# 读取基因列表和开始位点，长度为400bp
gene_file = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\candidate_list_2.txt"
gene_dict = dict()
with open(gene_file) as f:
    for line in f.readlines():
        l = line.split()
        gene_dict[l[0]] = l[1]

path_input = r"D:\python\3.MLST_SC\5.MLST\1.pangenome_analysis\2_genes_align"
path_output = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\1.candidate_fragment"
gene_allels =pd.DataFrame()
for gene, start_site in gene_dict.items():
    file_path = os.path.join(path_input, f"{gene}.aln.fas")
    s = int(start_site)
    df_alin = align_df(file_path).iloc[:, s: s+400]
    ds_allel = allel(df_alin)
    gene_allels[gene] = ds_allel
    seq = Seq(consensus_seq(df_alin))
    record = SeqRecord(id=gene, seq=seq, description="")
    file_output = os.path.join(path_output, f"{gene}.fa")
    SeqIO.write(record, file_output, "fasta")

gene_allels.to_excel("gene_allels.xlsx")