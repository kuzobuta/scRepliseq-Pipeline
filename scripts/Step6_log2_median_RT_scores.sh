#!/bin/bash

#scRepli-seq pipeline
#Step6 Computation of RT score at 200-kb window , 40-kb sliding

if [ $# -le 4 ] ; then
    echo "Usage: bash Step6_log2_median_RT_scores.sh [fragment_file] [out_dir] [ref_Rdata] [blacklist] [genome_file]"
    echo ""
    echo "Compute the log2 median profile for 200-kb window at 40-kb sliding interval"
    echo ""
    echo "arguments:"
    echo ""
    echo "fragment_file		fragment Rdata file (From Step3)"
    echo "out_dir			output directory ex) /Aneu_analysis/G1_control"
    echo "ref_Rdata			Merged control-G1 fragment Rdata file (From Step5_2)"
    echo "blacklist			black list file ex) ENCFF000KJP_name2id.bed etc.."
    echo "genome_file		genome file ex) hg19_female.fa.fai, mm9.chrom.size etc.."
    echo ""
    exit 0
fi

echo "#scRepli-seq Step6 computation of log2 median start: `date`"

rscript="/usr/local/bin/util/Step6_R_log2_median_RT_scores.R"
fragment=$1
outdir=$2
ref_Rdata=$3 
blacklist=$4
genome_file=$5

filename=`basename $1`
echo "#file name:${filename}"

Rscript --vanilla $rscript ${fragment} ${outdir} ${ref_Rdata} ${blacklist} ${genome_file}

echo "#scRepli-seq Step6 computation of log2 median start: `date`"
