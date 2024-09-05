from Bio import SeqIO

path_input = r"D:\python\3.MLST_SC\5.MLST\4.final_MLST\fragemet_alignments\atpB_2_fragment.aln.fas"
path_output = r"D:\python\3.MLST_SC\5.MLST\4.final_MLST\fragemet_alignments\atpB_2_fragements"
# record = SeqIO.read(path_file, "fasta")
record = SeqIO.parse(path_input, "fasta")

records = list(record)
for record in records:
    file_output = f"{path_output}/{record.id}.fa"
    SeqIO.write(record, file_output, "fasta")
