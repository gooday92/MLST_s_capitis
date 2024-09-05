import os


fa_list = []
gff_list = []


file_fa = open("fa.genomes.txt")
for line in file_fa.readlines():
    stain = line[0:-1].split(".")[0]
    fa_list.append(stain)
file_fa.close()

file_gff = open("gff.genomes.txt")
for line in file_gff.readlines():
    stain = line[0:-1].split(".")[0]
    gff_list.append(stain)
file_gff.close()

strain_list = []
for a in fa_list:
    if a not in gff_list:
        strain_list.append(a)
file = open("prokka.sh", "w")
for strain in strain_list:
    command = f" "
    file.write(command)
file.close()
