#!/bin/bash


PROJECT=../data
NAME=../data/fastq/hash1_1_Jvpna33.fq.gz

mkdir -pv $PROJECT/libs
NEWNAME=${NAME##*/}
NEWNAME=${NEWNAME%.fastq.gz}_trimmed.fastq.gz

#~/bin/bbmap/bbduk.sh -Xmx1g in=$NAME ref=~/bin/bbmap/resources/adapters.fa t=20 out=$PROJECT/libs/${NEWNAME}\
#		     ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=10 ftl=12
#fastqc $PROJECT/libs/${NEWNAME} -o $PROJECT/libs/

mkdir -pv $PROJECT/rna_align
DIR=$PROJECT/rna_align

~/bin/bbmap/bbmap.sh in=$PROJECT/libs/${NEWNAME} trimreaddescription=t  t=20 ref=$PROJECT/reference_sequences/FQ312003_wplasmids.fa\
		     k=8 ambig=best outm=$DIR/$NAME.sam 

