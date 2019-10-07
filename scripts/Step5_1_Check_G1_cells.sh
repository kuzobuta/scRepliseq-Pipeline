#!/bin/bash

#scRepli-seq pipeline
#Step5_1 Check the karyotype of the G1 cells by 6HMM in 500-kb bins

if [ $# -le 1 ] ; then
    echo "Usage: bash Step5_1_Check_G1_cells.sh [bin_file] [out_dir]"
    echo ""
    echo "6HMM calling by AneuFinder in 500-kb bins for checking the karyotype of G1 cells"
    echo ""
    echo "arguments:"
    echo ""
    echo "bin_file			bin data file (From Step3)"
    echo "out_dir			output directory ex) /Aneu_analysis/G1_control"
    echo ""
    exit 0
fi

echo "#scRepli-seq Step5_1 check G1 karyotype start: `date`"
filename=`basename $1`
echo "#file name:${filename}"

rscript="/usr/local/bin/util/Step5_R_G1_karyotype.R"
binfile=$1
out_dir=$2

Rscript --vanilla $rscript ${binfile} ${out_dir}

echo "#scRepli-seq Step5_1 check G1 karyotype start: `date`"
