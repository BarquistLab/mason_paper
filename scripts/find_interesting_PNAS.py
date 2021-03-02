import re
import json
import subprocess
from Bio import SeqIO
import os


# firstly get all essential locus tags in a list:
with open('../data/ESSENTIAL_GENES_salmonella.json', 'r') as f:
    ESSENTIAL_GENES = json.load(f)
essgenes = []
for i in ESSENTIAL_GENES['s_typhi']:
    egene = list(i.values())[1]
    essgenes += [egene]

# create new pdf with only ess. genes:

length_region = 10

# create gff with only essential salmonella genes:
with open("../data/FQ312003.1.gff") as origin_file:
    f = open("../data/essentialgenes.gff", "w")
    headers = re.compile(r"##")
    lt = re.compile(r".*;locus_tag=(SL1344_(P\d_)?\d{4})")

    for line in origin_file:
        if headers.match(line):
            f.write(line)
        elif lt.search(line) and "\tCDS\t" in line:
            locus_tag = lt.search(line)[1]
            if locus_tag in essgenes:
                line = line.split("\t")
                s = int(line[3])
                e = int(line[4])
                if line[6] == "+":
                    line[3] = str(s - (length_region - 3))
                    line[4] = str(s + length_region)
                else:
                    line[3] = str(e - length_region)
                    line[4] = str(e + (length_region - 3))
                f.write("\t".join(line))
    f.close()

# generate fasta of all essential gene start regions:
subprocess.call(["gffread", "-w", "../data/egenes.fa", "-g",
                 "../data/FQ312003_wplasmids.fa", "../data/essentialgenes.gff", "-F"])

# now change the fasta to locus tags headers:
subprocess.call("(sed -E 's/gene-([^ ]+).*/\\1/' ../data/egenes.fa > ../data/egenes1.fa)", shell=True)

# get 0 bp mismatches of all possible PNAs:
for seq_record in SeqIO.parse("../data/egenes1.fa", "fasta"):
    for i in range(9):
        newpna = SeqIO.SeqRecord(seq_record.seq[i:i+10], id=seq_record.id, description=str(i+1))
        savepath = "../data/pna_tests/" + seq_record.id + "_" + str(i+1)
        SeqIO.write(newpna, savepath + ".fasta", "fasta")
        subprocess.run(["seqmap", "0", savepath + ".fasta", "../data/egenes1.fa", savepath + "_result.tab",
                        "/output_all_matches", "/available_memory:200"])
        num_lines = sum(1 for line in open(savepath + "_result.tab"))
        if num_lines < 3:
            os.remove(savepath + "_result.tab")
        os.remove(savepath + ".fasta")
        print(newpna.seq)

double_pnas = {}

directory = "../data/pna_tests"
for filename in os.listdir(directory):
    file_name = os.path.join(directory, filename)
    with open(file_name, 'r') as f:
        next(f)
        for line in f:
            l = line.split("\t")
            try:
                if l[0] + "\t" + l[1] not in double_pnas[l[2]]:
                    double_pnas[l[2]] += [l[0] + "\t" + l[1]]
            except KeyError:
                double_pnas[l[2]] = [l[0] + "\t" + l[1]]

    with open("../data/new.tab", 'w') as f:
        for nt in double_pnas.keys():
            if len(double_pnas[nt]) > 1:
                line1 = nt + "\t" + double_pnas[nt][0] + "\n"
                line2 = nt + "\t" + double_pnas[nt][1] + "\n\n"
                f.write(line1)
                f.write(line2)


    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))


