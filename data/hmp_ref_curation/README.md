# Microbiome database

Date: August 23, 2022

Author: Jakob Jung

Here, we download the reference assemblies of microbiome datasets from the Human Microbiome Project (HMP) (Accession: PRJNA28331) and extract for each genome the regions around the start codon, which are potential off-targets. 



## microbiome.fasta & microbiome.gff

These were downloaded from [NCBIs HMP](https://www.ncbi.nlm.nih.gov/assembly?LinkName=bioproject_assembly_all&from_uid=28331) with the query  as described in https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/ on 7 of July 2022. Next, I generated a fasta file:

```bash
for NAME in fasta_files/*.fna.gz ; do  gunzip $NAME ; done
touch HMP.fasta
for NAME in fasta_files/*.fna ; do  cat $NAME >> HMP.fasta ; done

for NAME in gff_files/*.gff.gz ; do  gunzip $NAME ; done
touch HMP.gff
for NAME in gff_files/*.gff ; do  cat $NAME >> HMP.gff ; done
```

I add a line of code showing the length of all fastas using bioawk:

```bash
bioawk -c fastx '{ print $name, length($seq) }' < "HMP.fasta"  > "genelengths.tsv"
```



After creating the files, we used some linux commands to extract the regions around the start codons of all HMP genes of all HMP organisms in total . And We want to get the gene name/locus tag  into the GFF file into the first column to  :

```bash
grep -Pv "^#" HMP.gff | grep -P "\tgene\t" |awk 'BEGIN {OFS="\t"; FS="\t"}; { if ($7 == "+") {$5=$4+15;  $4=$4-30; if($4<1) $4=1} else if ($7 =="-"){ $4=$5-15;  $5=$5+30}; if($3 != "region") print $0}'  | grep -E "locus_tag="  | sed  -E 's/(([^\t]*)\t[^\t]*\t)[^\t]*(.*locus_tag=([^;]*).*)/\1\2_\4\3/' > start_regs_new.gff
```



I added an R script to handle genes which are too far in the end of a transcript, so not the whole 5"UTR can be captured:

```bash
Rscript modify_gff_HMP.R
```



Now we need to extract the start-region sequences of all of the entries in the gff to end up with a fasta file of very small sequences:

```bash
bedtools getfasta  -s -fi HMP.fasta -bed start_regs_new_mod.gff -name | awk '/^>/{f=!d[$1];d[$1]=1}f' > start_regions.fasta
```

Now we can use them for HMP data in mASON













