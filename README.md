# RNA-seq analysis to Jung et al. 2022

- Name of project: Design and off-target prediction for antisense oligomers targeting
  bacterial mRNAs with the MASON webserver
- Paper link to biorxiv: https://www.biorxiv.org/content/10.1101/2022.05.24.492283v1
- Data analysis: Jakob Jung
- Experiments: Linda Popella, Phuong Thao Do, Barbara Plaschke
- Supervision: Lars Barquist, Jörg Vogel
- Start: 2021
- End: November 2022

## Introduction

This project directory contains the analysis of RNA-Seq results obtained after Peptide-PNA challenge of *Salmonella enterica* subsp. *enterica*  serovar Typhimurium str. SL1344 with different PNA sequences. The procedure of RNA-seq experiments can be found in the according methods section of the [manuscript](https://www.biorxiv.org/content/10.1101/2022.05.24.492283v1).  

 

## Directory structure

The project is divided in 3 main subdirectories:

- data: contains all raw, intermediate and final data of the process.  
- analysis: contains analysis files, such as figures, plots, tables, etc. 
- scripts: contains scripts to process data from the data directory.

Each directory as well as subdirectories should have their own README.md file with information on the scripts, data and analysis files. 



## Workflow

Here I describe the workflow, which can be followed to reproduce the RNA-Seq results & plots from the article. 



### 1. Prerequisites

For running the whole RNA-seq analysis, one needs following packages/tools/software:

- BBMap (v38.84)
- R along with several packages from Bioconductor/CRAN 
- Linux for commands & bash scripts
- featureCounts (v2.0.1) from Subread package

 

### 2. Mapping

All raw FastQ files need to be put into the folder [./data/fastq](data/fastq). These fastq files can be accessed on GEO with accession number [GSE199542](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199542). Details on samples and upset of the experiment can be found in the methods section of the manuscript and fastQC quality statistics of the raw reads can be found in [./data/libs](./data/libs) . 

 To run the mapping, run the bash script [./scripts/trimm_map_BB.sh](./scripts/trimm_map_BB.sh) . The script loops through the fastq-files, trims of adapters using BBDuk, maps against the reference Salmonella genome (reference fasta and gff files can be found in [./data/reference_sequences/](./data/reference_sequences/)) and counts the reads mapped to coding sequences and sRNAs using featureCounts.

Trimming, mapping and counting statistics are stored in the log file [./analysis/logfile.log](./analysis/logfile.log) . The directory [./data/rna_align](./data/rna_align) includes all bam-alignment files as well as the count table [./data/rna_align/counttable.txt](./data/rna_align/counttable.txt) . This count will be imported into R for Differential expression analysis.



### 3. Differential expression analysis

Some prerequisite files need to be added beforehand such as the files [./data/PHOPQ.tsv](./data/PHOPQ.tsv)  and [./data/link_lt_gn.tab](./data/link_lt_gn.tab) . How they were generated is described in the README.md files of respective directories.

To run the differential expression analysis, run the R markdown script [./scripts/mason_RNASEQ.Rmd](./scripts/mason_RNASEQ.Rmd) . This outputs all figures of the manuscript, which are saved as PDF and/or SVG files to the [./analysis](./analysis) directory. It might take up to 10 minutes to run this script on a low-memory laptop. 