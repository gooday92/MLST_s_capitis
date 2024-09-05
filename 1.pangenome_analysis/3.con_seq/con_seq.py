from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os

def consensus_biopython(pathdir_input, pathfile_output):
    consensus_records = []
    for f in os.listdir(pathdir_input):
        if "aln.fas" in f:
            gene_name = f.split(".")[0]
            file_path = os.path.join(pathdir_input, f)
            # 读取alignment文件
            alignment = AlignIO.read(file_path, "fasta")
            # 总结alignment
            summary_align = AlignInfo.SummaryInfo(alignment)
            # 计算consensus序列，会引入X，有问题
            consensus = summary_align.dumb_consensus()
            # 以fasta格式保存序列
            record = SeqRecord(id=gene_name, seq=consensus, description="")
            consensus_records.append(record)
    SeqIO.write(consensus_records, pathfile_output, "fasta")


def consensus_pandas(pathdir_input, pathfile_output):
    consensus_records = []
    for f in os.listdir(pathdir_input):
        if "aln.fas" in f:
            gene_name = f.split(".")[0]
            file_path = os.path.join(pathdir_input, f)
            seq = Seq(consensus_seq(file_path))
            record = SeqRecord(id=gene_name, seq=seq, description="")
            path_out = fr"D:\python\3.MLST_SC\第三次MLST\2.calculate\1.con_seq\consesensus_seq\{gene_name}.fa"
            SeqIO.write(record, path_out, "fasta")
            # consensus_records.append(record)
            # print(consensus_records)
    # SeqIO.write(consensus_records, pathfile_output, "fasta")

def consensus_seq(file):
    alignment = AlignIO.read(file, "fasta")
    df = pd.DataFrame(alignment)
    seq = []
    for label, col in df.iteritems():
        # 取每一列占比最大的那个碱基
        cha = col.value_counts().index[0]
        seq.append(cha)
    conses_seq = "".join(seq)
    return conses_seq



pathdir_input = r"D:\python\3.MLST_SC\第三次MLST\1.pangenome_analysis\2_genes_align"
pathfile_output = r"D:\python\3.MLST_SC\第三次MLST\2.calculate\1.con_seq\genes_consensus_seq_pds.fasta"
# consensus_biopython(path_input, path_output)
consensus_pandas(pathdir_input, pathfile_output)




