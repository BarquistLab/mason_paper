# MFE predictions

Here I calculated mfes for all of the start regions, later used for a supplementary table in the paper.

I start by creating a gff with only the start regions (-30 until +15) of the genes:

```bash
#salmonella
grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t|\tgene\t" ../reference_sequences/FQ312003.1.gff | awk -F'\t' 'BEGIN { OFS="\t" } {if ($7=="-") {$4=$5-15 ; $5=$5+30 } else {$5=$4+15; $4=$4-30} print $0}' | sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;locus_tag=([^;]+).*)/\1\4\3/' | sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;gene=([^;]+).*)/\1\2;\4\3/'|grep ";locus_tag=" > salmonella_startregs.gff

#upec
grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t|\tgene\t" ../../../pna_rnaseq_experiments/upec_transcriptomics_03_21/data/reference_sequences/ecoli536.gff3 | awk -F'\t' 'BEGIN { OFS="\t" } {if ($7=="-") {$4=$5-15 ; $5=$5+30 } else {$5=$4+15; $4=$4-30} print $0}' | sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;locus_tag=([^;]+).*)/\1\4\3/' | sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*;gene=([^;]+).*)/\1\2;\4\3/'|grep ";locus_tag=" > upec_startregs.gff
```

now I wanna get the fasta files for all these sequences:

```bash
# extract sequences:
bedtools getfasta -s -fi ../reference_sequences/FQ312003_wplasmids.fa -bed "salmonella_startregs.gff" -name+ -fo "salm_startregs.fasta"

# upec
bedtools getfasta -s -fi ../../../pna_rnaseq_experiments/upec_transcriptomics_03_21/data/reference_sequences/ecoli536.fasta -bed "upec_startregs.gff" -name+ -fo "upec_startregs.fasta"
```



get MFEs:

```
RNAfold salm_startregs.fasta |sed -E 's/^>([^;:]+).*/\1/' |  perl -pe 's/\n/\t/g; s/>//; s/\s+/\t/; s/\(-|\( -/-/; s/\)\t$/\n/' | sed -E 's/([^\t]+)\t[AUGC\t\.\(\)  ]*(.+)/\1\t\2/' | sed 's/(  //' | uniq  > salm_mfe.tsv

#upec
RNAfold upec_startregs.fasta |sed -E 's/^>([^;:]+).*/\1/' |  perl -pe 's/\n/\t/g; s/>//; s/\s+/\t/; s/\(-|\( -/-/; s/\)\t$/\n/' | sed -E 's/([^\t]+)\t[AUGC\t\.\(\)  ]*(.+)/\1\t\2/' | sed 's/(  //' | uniq  > upec_mfe.tsv

#remove plots
rm -rf *.ps
```

