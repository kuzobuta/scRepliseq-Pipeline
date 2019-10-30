#!/bin/bash

#scRepli-seq pipeline
#Predict the SEQXE adaptor sequence from fastq file
util_dir="/usr/local/bin/util"

in_fastq=$1
out_dir=$2

echo "#scRepli-seq Predict the SEQXE adaptor sequence start: `date`"

mkdir -p ${out_dir}

fastq=`basename ${in_fastq}`
prefix=${out_dir}/${fastq%.fastq.gz}

#Extract first 50bps
seqtk trimfq -L 50 ${in_fastq} | gzip -c > ${prefix}.fastq.L50.gz

#Convert fastq into tab txt format
seqtk seq -A ${prefix}.fastq.L50.gz | perl ${util_dir}/fasta2tab.pl - | cut -f2 | sort | uniq -c | sort -nrk1 - > ${prefix}.tab.L50.sort.txt

#Grep the reads including 5'-CTGAAG-3' site & convert tab txt to fastq
grep CTGAAG ${prefix}.tab.L50.sort.txt | perl -ne 's/[ \t]+/\t/g;s/^\t//g; print' | perl ${util_dir}/tab2fastq.pl - | gzip -c - > ${prefix}.tab.L50.SEQXE.fastq.gz

#Check it by fastqc
fastqc ${prefix}.tab.L50.SEQXE.fastq.gz

echo "#scRepli-seq Predict the SEQXE adaptor sequence end: `date`"
echo "#Please check the FASTQC results and estimate the predicted SEQXE adaptor sequence"

