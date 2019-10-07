#!/bin/bash

#scRepli-seq pipeline
#Step7 Binarization of the data

if [ $# -le 4 ] ; then
    echo "Usage: bash Step7_Binarization.sh [binfile] [out_dir] [ref_Rdata] [genome_file] [binsize] [somy]"
    echo ""
    echo "Binarization of the bin data"
    echo ""
    echo "arguments:"
    echo ""
    echo "binfile		bin data file (From Step3)"
    echo "out_dir		output directory ex) /Aneu_analysis/HMM"
    echo "ref_Rdata		Merged control-G1 fragment Rdata file (From Step5_2)"
    echo "genome_file	genome file ex) hg19_female.fa.fai, mm9.chrom.size etc.."
    echo "binsize		binsize for the binarization ex) 80000, 100000"
    echo "somy			Mode of somy in HMM ex) 1-somy, 2-somy"
    echo ""
    echo "For Mid-S samples, you do not need to set the somy."
    echo ""
    exit 0
fi

echo "#scRepli-seq Step7 binarization start: `date`"

rscript="/usr/local/bin/util/Step7_R_Binarization.R"
binfile=$1
out_dir=$2
Ref_file=$3 #fragment file 
genome_file=$4
binsize=$5 #normally 80000 (80kb)
somy=$6 #For mid-S sample, you do not need to set this. Default: 2-somy

filename=`basename $1`
echo "#file name:${filename}"

Rscript --vanilla $rscript ${binfile} ${out_dir} ${Ref_file} ${genome_file} ${binsize} ${somy}

echo "#scRepli-seq Step7 binarization end: `date`"
