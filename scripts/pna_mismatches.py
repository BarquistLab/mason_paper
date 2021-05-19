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

# create new gff with only ess. genes:
length_region = 20

# create gff with only essential salmonella genes:
with open("../data/reference_sequences/FQ312003.1.gff") as origin_file:
    f_negenes = open("../data/mismatches_02_2021/nonessentialgenes.gff", "w")
    f_allgenes = open("../data/mismatches_02_2021/all_genes_startregions.gff", "w")

    headers = re.compile(r"##")
    lt = re.compile(r".*;locus_tag=(SL1344_(P\d_)?\d{4})")

    for line in origin_file:
        if headers.match(line):
            f_negenes.write(line)
        # write new start_region_gffs one for nonessential genes and one for all genes:
        elif lt.search(line) and "\tCDS\t" in line:
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
            f_allgenes.write("\t".join(line))
            if locus_tag not in essgenes:
                f_negenes.write("\t".join(line))

    f.close()

# generate fasta of all nonessential gene start regions:
subprocess.call(["gffread", "-w", "../data/mismatches_02_2021/nonessentialgenes.fa", "-g",
                 "../data/FQ312003_wplasmids.fa", "../data/mismatches_02_2021/nonessentialgenes.gff", "-F"])

# now change the fasta to locus tags headers:
subprocess.call("(sed -E 's/gene-([^ ]+).*/\\1/' ../data/mismatches_02_2021/nonessentialgenes.fa > "
                "../data/mismatches_02_2021/noness_genes_startregions.fa)", shell=True)

# generate fasta of all gene start regions:
subprocess.call(["gffread", "-w", "../data/mismatches_02_2021/all_genes_startregions.fa", "-g",
                 "../data/FQ312003_wplasmids.fa", "../data/mismatches_02_2021/all_genes_startregions.gff", "-F"])
# now change the fasta to locus tags headers:
subprocess.call("(sed -E 's/gene-([^ ]+).*/\\1/' ../data/mismatches_02_2021/all_genes_startregions.fa > "
                "../data/mismatches_02_2021/all_genes_startregion.fa)", shell=True)



# delete unmodified file:
os.remove("../data/mismatches_02_2021/nonessentialgenes.fa")
os.remove("../data/mismatches_02_2021/all_genes_startregions.fa")

# get 0 bp mismatches of all possible PNAs:
os.remove("../data/mismatches_02_2021/all_poss_pnas.fasta")
for seq_record in SeqIO.parse("../data/mismatches_02_2021/noness_genes_startregions.fa", "fasta"):
    for i in range(11, 17):
        newpna = SeqIO.SeqRecord(seq_record.seq[i:i+10], id=seq_record.id, description=str(i+1))
        newpna.id = seq_record.id + "_" + str(i+1)
        if GC(newpna.seq) > 50:
            with open("../data/mismatches_02_2021/all_poss_pnas.fasta", "a") as handle:
                SeqIO.write(newpna, handle, format="fasta")
        print(newpna.seq)


subprocess.run(
            ["seqmap", "2", "../data/mismatches_02_2021/all_poss_pnas.fasta",
             "../data/mismatches_02_2021/noness_genes_startregions.fa",
             "../data/mismatches_02_2021/all_pnas_2mm_result.tab", "/output_all_matches", "/available_memory:200"])

df_matches_raw = pd.read_csv("../data/mismatches_02_2021/all_pnas_2mm_result.tab", sep='\t')
df_matches_raw = df_matches_raw.sort_values("probe_id")
df_matches_raw.to_csv("../data/mismatches_02_2021/all_pnas_2mm_result.tab", sep="\t", index=False)

df_highexp = pd.read_csv("../data/strong_control.csv")
# get all locustags with highly expressed genes:
lt_highexp = df_highexp.iloc[:, 0].values

#  find pnas with low ots:
# find only 0,1 mm pairings:
df_matches = df_matches_raw[df_matches_raw["num_mismatch"].isin([0, 1])]
# remove duplicates to get only no-mismatch regions:
df_matches = df_matches.drop_duplicates(subset=["trans_id"])
# use only highly expressed genes:
df_matches = df_matches[df_matches["trans_id"].isin(lt_highexp)]
# remove mismatches in ess genes
df_matches = df_matches[~df_matches["trans_id"].isin(essgenes)]
# find only 0 mm pairings:
df_matches = df_matches[df_matches["num_mismatch"] == 0]
df_matches
df_matches.to_csv("../data/mismatches_02_2021/all_no_1_mismatch_pnas.tab", sep="\t")

#  find pnas with many ots:
# find only 0,1 mm pairings:
df_matches = df_matches_raw[df_matches_raw["num_mismatch"].isin([0, 1])]
df_matches.to_csv("../data/mismatches_02_2021/1_mismatch_pnas.tab", sep="\t")
df_muchofft = pd.DataFrame(columns=df_matches.columns)
c = 0
for index, row in df_matches.iterrows():

    if index == 0:
        continue
    else:
        if df_matches.iloc[index-1, 3] == df_matches.iloc[index, 3]:

            c += 1
        else:
            print(c)
            if c < 3:
                print(row)
                print(index)
                df_muchofft = df_muchofft.append([df_matches.iloc[index-1, :]])
            c = 0

c = 0
one_zero_mm = df_matches.iloc[:, 3].tolist()
fewofft = []
for lt in range(len(one_zero_mm)):
    if lt == 0:
        continue
    else:
        if one_zero_mm[lt] == one_zero_mm[lt-1]:
            c += 1
        else:
            if 3 < c < 6:
                fewofft += [one_zero_mm[lt - 1]]
                print(one_zero_mm[lt - 1], c)
            c = 0

for i in fewofft:
    x = re.search("(SL1344.*)_.*", i)[1]
    if x in lt_highexp:

        print(i, len(df_matches_raw.loc[(df_matches_raw["num_mismatch"] == 2) & (df_matches_raw["probe_id"] == i)]))








df_muchofft.to_csv("../data/mismatches_02_2021/all_high_offt_pnas.tab", sep="\t")


# check single pna:
id = "SL1344_1163_17 17"
mismatches = 1
print(df_matches_raw.loc[(df_matches_raw["num_mismatch"] == mismatches) & (df_matches_raw["probe_id"] == id)])



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

