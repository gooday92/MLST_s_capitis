import os

dir_path = r"D:\python\3.MLST_SC\5.MLST\4.final_MLST\fragemet_alignments"
dir_out = r"D:\python\3.MLST_SC\5.MLST\4.final_MLST\fragemet_alignments2"
for f in os.listdir(dir_path):
    if f.endswith("fa"):
        gene = f.split(".")[0]
        path = os.path.join(dir_path, f)
        with open(path) as file:
            file_new = open(f"{gene}_aln.fas", "w")
            for line in file.readlines():
                if line.startswith(">"):
                    line_new = line.split(";")[0]
                    line_new = f"{line_new}\n"
                else:
                    line_new = line
                file_new.write(line_new)
            file_new.close()