import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import os
from icecream import ic
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def align_df(path):
    alignment = AlignIO.read(file_path, "fasta")
    genes_list = []
    for record in alignment:
        genes_list.append(record.id)
    df = pd.DataFrame(alignment, index=genes_list)
    return df



genes_list = dict()
with open("fragment_start.txt") as file:
    for line in file.readlines():
        strain = line.split()[0]
        start = line.split()[1]
        genes_list[strain] = start


align_dir = r"D:\python\3.MLST_SC\5.MLST\1.pangenome_analysis\2_genes_align"
output_dir = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\5.final_MLST\fa_alleles"
alleles_seq = dict()
for gene in genes_list:
    file_path = os.path.join(align_dir, f"{gene}.aln.fas")
    df = align_df(file_path)
    start = int(genes_list[gene])
    df_slide = df.iloc[:, start:start+400]
    alleles_seq[gene] = dict()
    for gene_, row in df_slide.iterrows():
        type_seq = "".join(row)
        gene_ = gene_.split(";")[0]
        alleles_seq[gene][gene_] = type_seq
df_alleles_seq = pd.DataFrame(alleles_seq)
# ic(df_alleles_seq)


df_alleles_types = pd.read_excel(r"strain_ST.xlsx", index_col=0)
df_alleles_types = df_alleles_types[genes_list]
df_alleles_seq_types = pd.concat([df_alleles_seq, df_alleles_types], axis=1)
# ic(df_alleles_seq_types)
# ic(df_alleles_seq_types.columns)

for gene in genes_list:
    a = df_alleles_seq_types[gene]
    ic(a)
    for strain, row in a.iterrows():
        seq_ = Seq(row.iloc[0])
        type = row.iloc[1]
        id = f"{gene}_{type}"
        file_path = os.path.join(output_dir, f"{id}.fa")
        record = SeqRecord(seq=seq_, id=id, description="")
        SeqIO.write(record, file_path, "fasta")


