#!/bin/bash

if [ $# -le 2 ] ; then
    echo ""
    echo "Usage: bash Fragment2bin.sh <fragment file directory> <out_dir> <genome file> <bin size> 

#<bin size> is a option
#If <bin size> is not specified, default binsizes are selected as 40,80,100,200,500kb & 1Mb
#Input fragment R file is XXX_mapq10_blacklist_fragment.Rdata"
    echo ""
    exit 0
fi

rscript="/usr/local/bin/util/fragment2bin.R"
frag_dir=$1
out_dir=$2
genome_file=$3
binsize=$4

Rscript $rscript ${frag_dir} ${out_dir} ${genome_file} ${binsize}
