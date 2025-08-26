#!/bin/bash
#scRepli-seq pipeline
#Step2 Mapping & Mark duplicates
#v250825

if [ $# -le 4 ] ; then
    echo "Usage: bash Step2_Mapping.sh [Fastq_file] [bwa_index_file] [genome_name] [thread] [out_bam_file_name]"
    echo ""
    echo "Mapping, Marking the duplicated reads by picard following CleanSam, and Sorting by samtools"
    echo ""
    echo "arguments:"
    echo ""
    echo "Fastq_file		Input fastq file"
    echo "bwa_index_file	Index file for bwa mapping"
    echo "genome_name		genome short name ex) mm9, hg19, mm9_female etc.."
    echo "Thread			The number of threads for bwa mapping"
    echo "out_bam_file_name	Output bam file name"
    echo ""
    echo "Processed bam file name will be"
    echo "[out_bam_file_name].[genome_name].clean_srt_markdup.bam"
    echo ""
    exit 0
fi

echo "#scRepli-seq Step2 Mapping start: $(date)"


#Step2.1 Mapping
FASTQ=$1
index=$2
genome=$3 #Ex: mm9, hg19
THREAD=$4 #It depends on your machine
OUTNAME=$5

echo "#file name:${filename}"

if [ -z ${index}.bwt ]
then
echo "# Bwa index was not found. Please generate the index.."
exit 0
else
echo "# BWA index file was found. Skip building the index files.."
fi

bwa aln -t ${THREAD} ${index} ${FASTQ} | bwa samse ${index} - ${FASTQ} | samtools view -Sb - > ${OUTNAME}.${genome}.bam

echo "#scRepli-seq Step2 Mapping end: `date`"
echo "#scRepli-seq Step2 Clean & Sort start: `date`"


#Step2.2 Clean & Sort & Markduplicates

input_bam=${OUTNAME}.${genome}.bam
out_bam1=${OUTNAME}.${genome}.clean.bam
out_bam1_srt=${OUTNAME}.${genome}.clean_srt.bam
out_bam1_srt_markdup=${OUTNAME}.${genome}.clean_srt_markdup.bam
out_bam1_srt_markdup_met=${OUTNAME}.${genome}.clean_srt_markdup_met.txt

picard CleanSam I=$input_bam O=$out_bam1

#Sort the bam file
samtools sort $out_bam1 > $out_bam1_srt
samtools index $out_bam1_srt

echo "#scRepli-seq Step2 Clean & Sort end: `date`"
echo "#scRepli-seq Step2 Markduplicates start: `date`"


#Mark the duplicate reads & index
picard MarkDuplicates I=$out_bam1_srt O=$out_bam1_srt_markdup \
METRICS_FILE=$out_bam1_srt_markdup_met \
REMOVE_DUPLICATES=false

samtools index $out_bam1_srt_markdup

#Remove the redundant bam files
rm $out_bam1 $out_bam1_srt ${out_bam1_srt}.bai

echo "#scRepli-seq Step2 Markduplicates end: `date`"

