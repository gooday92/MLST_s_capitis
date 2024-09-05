from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from icecream import ic

def phy2fa(input_file, output_file):
    alignment = []
    with open(input_file) as file:
        for line in file.readlines():
            if len(line) > 1:
                strain = line.split()[0]
                seq = line.split()[1]
                if len(seq) > 100:
                    record = SeqRecord(seq=Seq(seq), id=strain, description="")
                    ic(record)
                    alignment.append(record)
    alignment = MultipleSeqAlignment(alignment)
    AlignIO.write(alignment, output_file, "fasta")

file_path = r"D:\python\3.MLST_SC\5.第五次MLST\2.groups\fastbaps\core_gene_alignment.aln.varsites.phy"
out = r"all_core_genes.fa"
phy2fa(file_path, out)
