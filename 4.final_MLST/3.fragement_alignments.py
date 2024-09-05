import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
import os
from icecream import ic
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def main(dir_input, file, dir_out):
    # 解析文件获得基因列表和起始位点
    gene_dict = dict()
    with open(file) as f:
        for line in f.readlines():
            gene = line.split()[0]
            start = line.split()[1]
            gene_dict[gene] = start
    # 读取相应的文件截取相应的alignment片段
    for gene,start in gene_dict.items():
        start = int(start)
        # ic(gene, start)
        file_input = f"{dir_input}/{gene}.aln.fas"
        alignment = AlignIO.read(file_input, "fasta")
        alignment_new = alignment[:, start:start+400]
        file_new = fr"{dir_out}/{gene}_fragment.aln.fa"
        AlignIO.write(alignment_new, file_new, "fasta")


dir_input = r"D:\python\3.MLST_SC\5.MLST\1.pangenome_analysis\2_genes_align"
file_input = r"D:\python\3.MLST_SC\5.MLST\4.final_MLST\fragment_start.txt"
dir_out = r"D:\python\3.MLST_SC\5.MLST\4.final_MLST\fragemet_alignments"
main(dir_input, file_input, dir_out)
