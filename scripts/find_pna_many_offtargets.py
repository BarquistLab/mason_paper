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
with open("../data/reference_sequences/FQ312003.1.gff") as origin_file:
    f_ess = open("../data/essentialgenes_start.gff", "w")
    f_noness = open("../data/nonessentialgenes_start.gff", "w")
    headers = re.compile(r"##")
    lt = re.compile(r".*;locus_tag=(SL1344_(P\d_)?\d{4})")
    gene = re.compile(r".*;gene=([^;]+)")

    for line in origin_file:
        if headers.match(line):
            f_noness.write(line)
            continue
        elif lt.search(line) and "\tCDS\t" in line:
            if gene.search(line):
                genename = gene.search(line)[1]
            else:
                genename = ""

            locus_tag = lt.search(line)[1]
            line = line.split("\t")
            s = int(line[3])
            e = int(line[4])
            if line[6] == "+":
                line[3] = str(s - (length_region - 3))
                line[4] = str(s + length_region)
            else:
                line[3] = str(e - length_region)
                line[4] = str(e + (length_region - 3))
            if locus_tag in essgenes:
                f_ess.write("\t".join(line))
            f_noness.write("\t".join(line))
    f_ess.close()
    f_noness.close()

# generate fasta of all essential gene start regions:
subprocess.call(["gffread", "-w", "../data/essentialgenes_start.fa", "-g",
                 "../data/FQ312003_wplasmids.fa", "../data/essentialgenes_start.gff", "-F"])

# generate fasta of all noness' start regions:
subprocess.call(["gffread", "-w", "../data/nonessentialgenes_start.fa", "-g",
                 "../data/FQ312003_wplasmids.fa", "../data/nonessentialgenes_start.gff", "-F"])

# now change the fasta to locus tags headers:
subprocess.call("(sed -E 's/gene-([^ ]+).*/\\1/' ../data/essentialgenes_start.fa > ../data/essgenes_start.fa)",
                shell=True)
subprocess.call("(sed -E 's/gene-([^ ]+).*/\\1/' ../data/nonessentialgenes_start.fa > "
                "../data/nonessgenes_start.fa)", shell=True)

# get 0 bp mismatches of all possible PNAs:
for seq_record in SeqIO.parse("../data/essgenes_start.fa", "fasta"):
    for i in range(9):
        newpna = SeqIO.SeqRecord(seq_record.seq[i:i+10], id=seq_record.id, description=str(i+1))
        savepath = "../data/pna_mismatches/" + seq_record.id + "_" + str(i+1)

        # add

        SeqIO.write(newpna, savepath + ".fasta", "fasta")
        subprocess.run(["seqmap", "0", savepath + ".fasta", "../data/nonessgenes_start.fa", savepath + "_result.tab",
                        "/output_all_matches", "/available_memory:200", "/forward_strand"])
        num_lines = sum(1 for line in open(savepath + "_result.tab"))
        if num_lines < 3:
            os.remove(savepath + "_result.tab")
        os.remove(savepath + ".fasta")
        print(newpna.seq)

double_pnas = {}

directory = "../data/pna_mismatches"
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

with open("../data/pna_mismatches/pna_mm.tab", 'w') as f:
    for nt in double_pnas.keys():
        if len(double_pnas[nt]) > 1:
            for gene in double_pnas[nt]:
                line1 = nt + "\t" + gene + "\n"
                f.write(line1)
            f.write("\n")
