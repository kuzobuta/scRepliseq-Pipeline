#!/bin/bash

#scRepli-seq pipeline
#Step4 Check the Quality of the data in 200-kb bins by MAD scoring
#Compute the MAD scores for all samples in "bins_data_dir" directory

if [ $# -le 2 ] ; then
    echo "Usage: bash Step4_MAD_score_QC.sh [bins_data_dir] [out_dir] [out_name]"
    echo ""
    echo "MAD score calculation in 200-kb bins and export into the file"
    echo ""
    echo "arguments:"
    echo ""
    echo "bins_data_dir		bins data directory (From Step3)"
    echo "out_dir			output directory"
    echo "out_name			output name ([out_name]_MAD_scores_log2.txt)"
    echo ""
    echo "Output file includes the MAD socres for individual file"
    echo "Rdata_file_Name	MAD_scores"
    echo ""
    exit 0
fi

echo "#scRepli-seq Step4 MAD scoring start: `date`"

rscript="/usr/local/bin/util/Step4_R_MAD_score.R"
bins_data_dir=$1
out_dir=$2
out_name=$3

Rscript --vanilla $rscript ${bins_data_dir} ${out_dir} ${out_name}

echo "#scRepli-seq Step4 MAD scoring end: `date`"


