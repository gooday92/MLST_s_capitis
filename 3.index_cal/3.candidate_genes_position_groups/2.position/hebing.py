from Bio import SeqIO
import os


path_input = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\1.candidate_fragment"
path_output = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\2.position"
file_name = "candidate_genes.fa"

records = []
for file in os.listdir(path_input):
    file_fa = os.path.join(path_input,file)
    record = SeqIO.read(file_fa,"fasta")
    records.append(record)
file_output = os.path.join(path_output,file_name)
SeqIO.write(records, file_output, "fasta")