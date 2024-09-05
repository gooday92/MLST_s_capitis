import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


path_input =rf"D:\python\3.MLST_SC\第三次MLST\2.calculate\1.con_seq\consesensus_seq"
path_output = rf"D:\python\3.MLST_SC\第三次MLST\2.calculate\1.con_seq\consesensus_seq_no_"


for file in os.listdir(path_input):
    file_path = os.path.join(path_input, file)
    file_output = os.path.join(path_output, file)
    record = SeqIO.read(file_path, "fasta")
    seq_new = Seq(str(record.seq).replace("-", ""))
    record_new = SeqRecord(seq_new, id=record.id, description="")
    SeqIO.write(record_new, file_output, "fasta")


