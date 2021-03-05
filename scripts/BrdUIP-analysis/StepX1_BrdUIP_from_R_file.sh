#!/bin/bash

#scRepli-seq pipeline for BrdU-IP analysis
#From R fragment file generated by Step3_load_mapped_reads.sh 
#This script will generate bedgraph file of BrdUIP 
#ver. 210304

if [ $# -le 4 ] ; then
    echo "Usage: bash StepX1_BrdUIP_from_bam_file.sh [Early_S] [Late_S] [name] [out_dir] [blacklist] [genome_file]"
    echo ""
    echo "Load the mapped reads, convert into Rdata format, and binning the data in Rdata format"
    echo ""
    echo "arguments:"
    echo ""
    echo "Early_S			Input fragment R file (Early-S fraction generated by Step3_load_mapped_reads.sh)"
    echo "Late_S			Input fragment R file (Late-S fraction generated by Step3_load_mapped_reads.sh)"
    echo "name				Output file name"
    echo "out_dir			Output directory (generating bins & fragment folder)"
    echo "blacklist			black list file ex) ENCFF000KJP_name2id.bed etc.."
    echo "genome_file		genome file ex) hg19_female.fa.fai, mm9.chrom.size etc.."
    echo ""
    echo "Output file name will be"
    echo "[name]_IP_Aneu_rmdup_w200ks40k_Percent_q0.05.bedGraph"
    echo ""
    exit 0
fi

echo "#Repli-seq StepX_IP analysis of BrdUIP from R file data start: `date`"

rscript="/usr/local/bin/util/repliseq.r"
EarlyS=$1
LateS=$2
name=$3 
out_dir=$4
blacklist=$5
genome_file=$6

filename=`basename $bamfile`
echo "#file name: ${EarlyS} and ${LateS}"

Rscript --vanilla $rscript ${EarlyS} ${LateS} ${name} ${out_dir} ${blacklist} ${genome_file}

echo "#Repli-seq StepX_IP analysis of BrdUIP data from R file end: `date`"
