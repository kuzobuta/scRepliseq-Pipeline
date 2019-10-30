#!/bin/bash

#scRepli-seq pipeline
#Step1.1 Trim illumina adaptor sequence

if [ $# -le 2 ] ; then
    echo "Usage: bash Step1_Trim_fastq.sh [index_seq] [in_fastq] [out_fastq1] [SEQXE] [out_fastq2]"
    echo ""
    echo "Trimming the fastq file by cutadapt program using index adaptor & SEQXE primer sequence"
    echo ""
    echo "arguments:"
    echo ""
    echo "index_seq     Illumina index adaptor sequence"
    echo "in_fastq		The path of fastq file"
    echo "out_fastq1	The path of adaptor trimmed fastq file"
    echo "SEQXE			SEQXE adaptor sequence (optional)"
    echo "out_fastq2	The path of adaptor & SEQXE-primer trimmed fastq file (optional)"
    echo ""
    echo "For TruSeq Single Indexes, you can use index_seq as "
    echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
    echo ""
    exit 0
fi

echo "#scRepli-seq Step1 Adaptor trimming	start: `date`"

#This index sequence is for Illumina TruSeq Single Indexes
index_seq=$1
#index_seq="GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

in_fastq=$2
out_fastq1=$3 #ex: QUERY.adapter_filtered.fastq.gz 

filename=`basename $in_fastq`
echo "#file name:${filename}"


cutadapt -b ${index_seq} -o ${out_fastq1} ${in_fastq}

echo "#scRepli-seq Step1 Adaptor Trimming	end: `date`"
echo "#scRepli-seq Step1 SEQXE Trimming		start: `date`"


#Step1.2 Trim SEQXE primer
SEQXE=$4
out_fastq2=$5 #ex: QUERY.adapter_filtered2.fastq.gz 

if [ -z "$SEQXE" ]
then

echo "#SEQXE Primer trimming is not performed.."

else

min_length=30
cutadapt -b ${SEQXE} -e 0.09 -O 19 -m ${min_length} -o ${out_fastq2} ${out_fastq1}

fi

echo "#scRepli-seq Step1 Whole trimming steps	end: `date`"
