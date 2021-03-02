import re
import json
import subprocess
from Bio import SeqIO
import os
from Bio.SeqUtils import GC
import pandas as pd


# firstly get all essential locus tags in a list:
with open('../data/ESSENTIAL_GENES_salmonella.json', 'r') as f:
    ESSENTIAL_GENES = json.load(f)
essgenes = []
for i in ESSENTIAL_GENES['s_typhi']:
    egene = list(i.values())[1]
    essgenes += [egene]

# create new pdf with only ess. genes:

length_region = 20

# create gff with only essential salmonella genes:
with open("../data/FQ312003.1.gff") as origin_file:
    f = open("../data/mismatches_02_2021/nonessentialgenes.gff", "w")
    headers = re.compile(r"##")
    lt = re.compile(r".*;locus_tag=(SL1344_(P\d_)?\d{4})")

    for line in origin_file:
        if headers.match(line):
            f.write(line)
        elif lt.search(line) and "\tCDS\t" in line:
            locus_tag = lt.search(line)[1]
            if locus_tag not in essgenes:
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

# generate fasta of all nonessential gene start regions:
subprocess.call(["gffread", "-w", "../data/mismatches_02_2021/nonessgenes.fa", "-g",
                 "../data/FQ312003_wplasmids.fa", "../data/mismatches_02_2021/nonessentialgenes.gff", "-F"])

# now change the fasta to locus tags headers:
subprocess.call("(sed -E 's/gene-([^ ]+).*/\\1/' ../data/mismatches_02_2021/nonessgenes.fa > "
                "../data/mismatches_02_2021/negenes1.fa)", shell=True)

# get 0 bp mismatches of all possible PNAs:
for seq_record in SeqIO.parse("../data/mismatches_02_2021/negenes1.fa", "fasta"):
    for i in range(11, 14):

        newpna = SeqIO.SeqRecord(seq_record.seq[i:i+10], id=seq_record.id, description=str(i+1))
        newpna.id = seq_record.id + "_" + str(i+1)
        if GC(newpna.seq) > 50:
            with open("../data/mismatches_02_2021/all_poss_pnas.fasta", "a") as handle:
                SeqIO.write(newpna, handle, format="fasta")
        print(newpna.seq)


subprocess.run(
            ["seqmap", "0", "../data/mismatches_02_2021/all_poss_pnas.fasta", "../data/mismatches_02_2021/negenes1.fa",
             "../data/mismatches_02_2021/all_pnas_result.tab", "/output_all_matches", "/available_memory:200"])

df_matches = pd.read_csv("../data/mismatches_02_2021/all_pnas_result.tab", sep='\t')
df_highexp = pd.read_csv("../data/strong_control.csv")
# get all locustags with highly expressed genes:
lt_highexp = df_highexp.iloc[:500,0].values

# remove duplicates to get only no-mismatch regions:
df_matches = df_matches.drop_duplicates(subset=["trans_id"])
# use only highly expressed genes:
df_matches = df_matches[df_matches["trans_id"].isin(lt_highexp)]

df_matches.to_csv("../data/mismatches_02_2021/all_nomismatch_pnas.tab" , sep="\t")






# get 0 bp mismatches of all possible PNAs:
for seq_record in SeqIO.parse("../data/mismatches_02_2021/negenes1.fa", "fasta"):
    for i in range(11, 14):
        newpna = SeqIO.SeqRecord(seq_record.seq[i:i+10], id=seq_record.id, description=str(i+1))
        savepath = "../data/mismatches_02_2021/pna_tests/" + seq_record.id + "_" + str(i+1)
        SeqIO.write(newpna, savepath + ".fasta", "fasta")
        subprocess.run(["seqmap", "0", savepath + ".fasta", "../data/mismatches_02_2021/negenes1.fa", savepath + "_result.tab",
                        "/output_all_matches", "/available_memory:200"])
        num_lines = sum(1 for line in open(savepath + "_result.tab"))
        if num_lines < 6 or GC(newpna.seq) < 50:
            os.remove(savepath + "_result.tab")
        os.remove(savepath + ".fasta")
        print(newpna.seq)

double_pnas = {}

directory = "../data/mismatches_02_2021/pna_tests"
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

    with open("../data/mismatches_02_2021/new.tab", 'w') as f:
        for nt in double_pnas.keys():
            if len(double_pnas[nt]) > 1:
                line1 = nt + "\t" + double_pnas[nt][0] + "\n"
                line2 = nt + "\t" + double_pnas[nt][1] + "\n\n"
                f.write(line1)
                f.write(line2)


    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))

