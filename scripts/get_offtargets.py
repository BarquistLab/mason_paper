"""
Here we create gff and fasta file with only start regions to search for off-targets. We use seqmap to identify them.
"""

import re
import subprocess
import os
from Bio.Seq import Seq


# create new gff with only ess. genes:
length_region = 15

# create gff with only essential salmonella genes:
with open("../data/reference_sequences/FQ312003.1.gff") as origin_file:
    f_start_regions = open("../data/mismatches_02_2021/start_regions.gff", "w")

    headers = re.compile(r"##")
    mRNA = re.compile(r"\tCDS\t|\trRNA\t|\ttRNA\t|\tsRNA\t")

    for line in origin_file:
        if headers.match(line):
            f_start_regions.write(line)
        # write new start_region_gffs one for nonessential genes and one for all genes:
        elif mRNA.search(line):
            line = line.split("\t")
            s = int(line[3])
            e = int(line[4])
            if line[6] == "+":
                line[3] = str(s - 15)
                line[4] = str(s + 11)
            else:
                line[3] = str(e - 11)
                line[4] = str(e + 15)
            f_start_regions.write("\t".join(line))
    f_start_regions.close()

# generate fasta of all nonessential gene start regions:
subprocess.call(["gffread", "-w", "../data/mismatches_02_2021/start_regions1.fa", "-g",
                 "../data/reference_sequences/FQ312003_wplasmids.fa", "../data/mismatches_02_2021/start_regions.gff",
                 "-F"])

# now change the fasta to locus tags headers (also run in command line!):
subprocess.call("sed -E 's/gene-([^ ]+).*/\1/' ../data/mismatches_02_2021/start_regions1.fa | "
                "sed -E 's/.*(locus_tag|product)=([^ \n]+).*/>\2/' > "
                "../data/mismatches_02_2021/start_regions.fa", shell=True)



# delete unmodified file:
os.remove("../data/mismatches_02_2021/start_regions1.fa")

# now I manually assembled all pna target sequences using:
# my_pna = Seq("CTCATCTGTC")
# print(my_pna.reverse_complement())
# for each PNA to get the target sequence. it is saved as "pna_targets.fasta"

# I created a fasta file from out PNA sequences for seqmap:

with open("../data/mismatches_02_2021/PNA_sequences.csv", "r") as f:
    with open("../data/mismatches_02_2021/PNA_sequences.fasta", "w") as fasta:
        next(f)
        for line in f:
            line = line[:-1]
            s = line.split("\t")
            fasta.write(">" + s[0] + "_" + s[1] + "\n")
            sequence = Seq(s[2])
            s_revcomp = sequence.reverse_complement()
            fasta.write(str(s_revcomp) + "\n")
            print(s[1])
            print(s_revcomp)


# now I check off-targets with up to 2 mm in target regions ....
# (I ran it in terminal, because subprocess didnt work somehow...)
# ran same with 3 mismatches!

subprocess.run(
    "seqmap 2 ../data/mismatches_02_2021/PNA_sequences.fasta ../data/mismatches_02_2021/start_regions.fa "
    "../data/mismatches_02_2021/pna_2mm_startregions.tab "
    "/output_all_matches /available_memory:200 /forward_strand",
    shell=True)

