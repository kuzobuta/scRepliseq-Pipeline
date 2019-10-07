#!/bin/bash

#scRepli-seq pipeline
#Step5_2 Merge G1 cells

if [ $# -le 1 ] ; then
    echo "Usage: bash Step5_2_Merge_G1_cells.sh [Merged R data file name] [fragment R file1] [fragment R file2] .."
    echo ""
    echo "Merge the fragment R data file into single fragment R file"
    echo ""
    echo ""
    exit 0
fi

echo "#scRepli-seq Step5_2 Merge G1 cells start: `date`"

rscript="/usr/local/bin/util/Merge_fragment_Rdata.R"
args=($@)
files=(${args[@]:1})

echo "#Merged files are:"
echo "#Total ${#files[@]} files"
for file in ${files[@]}
do
file=`basename $file`
echo "# $file"
done

Rscript --vanilla $rscript ${files[@]} -o ${args[0]}

echo "#scRepli-seq Step5_2 Merge G1 cells end: `date`"

