#!/bin/bash
#scRepli-seq pipeline
#Step1.1 Trim illumina adaptor sequence
#Multi-thread mode
#v250825

if [ $# -le 4 ] ; then
    echo "Usage: bash Step1_Trim_fastq.sh [index_seq] [threads] [in_fastq] [outdir] [SEQXE]"
    echo ""
    echo "Trimming the fastq file by cutadapt program using index adaptor & SEQXE primer sequence"
    echo ""
    echo "arguments:"
    echo ""
    echo "index_seq     Illumina index adaptor sequence"
    echo "threads       CPUs for cutadapt"
    echo "in_fastq		The path of fastq file"
    echo "outdir    	Outdir path"
    echo "SEQXE			SEQXE adaptor sequence (optional)"
    echo ""
    echo "For TruSeq Single Indexes, you can use index_seq as "
    echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
    echo "For SEQXE, you may use the following seq "
    echo "TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN"
    echo ""
    exit 0
fi

echo "#scRepli-seq Step1 Adaptor trimming	start: `date`"

#This index sequence is for Illumina TruSeq Single Indexes
index_seq=$1
#index_seq="GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

# CPUs for cutadapt
threads=$2

in_fastq=$3
filename=`basename $in_fastq`
echo "#file name:${filename}"

outdir=$4
mkdir -p $outdir

prefix=${filename%.fastq.gz}
out_fastq1="${outdir}/${prefix}.adapter_filtered.fastq.gz"

if [-e "${out_fastq1}" ]; then

echo "#"
echo "# Skip adapter trimming.. Found adapter filtered fastq for ${filename}" 
echo "#"

else

echo "#scRepli-seq Step1 Adaptor trimming   start: $(date)"
cutadapt -j $threads -b ${index_seq} -o ${out_fastq1} ${in_fastq}
echo "#scRepli-seq Step1 Adaptor Trimming	  end: $(date)"

fi


#Step1.2 Trim SEQXE primer
SEQXE=$5
out_fastq2=$5 #ex: QUERY.adapter_filtered2.fastq.gz 

if [ -z "$SEQXE" ]
then

echo "#SEQXE Primer trimming is not performed.."

else

if [ -e "${out_fastq2}" ]; then
echo "#"
echo "# Skip SEQXE trimming.. Found SEQXE filtered fastq for ${filename}"
echo "#"

else

echo "#scRepli-seq Step1 SEQXE Trimming     start: $(date)"
min_length=30
cutadapt -b ${SEQXE} -e 0.09 -O 19 -m ${min_length} -o ${out_fastq2} ${out_fastq1}
echo "#scRepli-seq Step1 Whole trimming steps end: $(date)"

fi

fi
