#!/bin/bash
#scRepli-seq test run for checking scripts
#v260225
set -euo pipefail

echo "
# Start test run:  $(date)"

start=$(date +%s)

# Setup path
bash Setup_path_for_local.sh

# ---------------------
# 00. Setup directories 
# ---------------------

OUTDIR="test_run"
bash scripts/Step0_setup_directories.sh $OUTDIR

echo "# Done step0:      $(date)"

# ---------------------
# FastQC 
# ---------------------

fastqc -o "$OUTDIR/fastqc" sample_data/fastq/Sample_P285_09_1_R1_100k.fastq.gz &> "$OUTDIR/logs/fastqc.log"
echo "# Done FastQC:     $(date)"

# ---------------------
# 01. Trim fastq
# ---------------------

in_fastq="sample_data/fastq/Sample_P285_09_1_R1_100k.fastq.gz"
index_seq="GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
out_fastq1="$OUTDIR/trim_fastq/QUERY.adapter_filtered.fastq.gz"
SEQXE_PRIMER="TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN"
out_fastq2="$OUTDIR/trim_fastq/QUERY.adapter_filtered2.fastq.gz"

bash scripts/Step1_Trim_fastq.sh $index_seq $in_fastq $out_fastq1 $SEQXE_PRIMER $out_fastq2 &> "$OUTDIR/logs/Step1.log"

echo "# Done step1:      $(date)"

# ---------------------
# 02-1. Map
# ---------------------

mkdir -p $OUTDIR/references
wget -q -O $OUTDIR/references/chrM.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chrM.fa.gz 2>/dev/null
bwa index $OUTDIR/references/chrM.fa.gz &> $OUTDIR/references/bwa.index.log

out_fastq2="$OUTDIR/trim_fastq/QUERY.adapter_filtered2.fastq.gz"
THREAD=2 #It depends on your machine.
index="$OUTDIR/references/chrM.fa.gz"
FASTQ=${out_fastq2}
OUTNAME="$OUTDIR/bam/QUERY.adapter_filtered2"
genome_name="hg19_chrM" #short name of your reference genome (i.e. hg19, mm10, etc.)

bash scripts/Step2_Mapping.sh $out_fastq2 $index $genome_name $THREAD $OUTNAME &> "$OUTDIR/logs/Step2.log"

# ---------------------
# 02-2. clean bam & markdup
# ---------------------

input_bam="sample_data/bam/QUERY.adapter_filtered2.hg19_female.bam"
out_bam1="$OUTDIR/bam/QUERY.adapter_filtered2.hg19_female.clean.bam"
out_bam1_srt="$OUTDIR/bam/QUERY.adapter_filtered2.hg19_female.clean_srt.bam"
out_bam1_srt_markdup="$OUTDIR/bam/QUERY.adapter_filtered2.hg19_female.clean_srt_markdup.bam"
out_bam1_srt_markdup_met="$OUTDIR/bam/QUERY.adapter_filtered2.hg19_female.clean_srt_markdup.metrices.txt"
picard CleanSam I=$input_bam O=$out_bam1 &> "$OUTDIR/logs/picard.CleanSam.log"

samtools sort $out_bam1 > $out_bam1_srt
samtools index $out_bam1_srt

picard MarkDuplicates I=$out_bam1_srt O=$out_bam1_srt_markdup \
METRICS_FILE=$out_bam1_srt_markdup_met \
REMOVE_DUPLICATES=false \
&> "$OUTDIR/logs/picard.MarkDup.log"

samtools index $out_bam1_srt_markdup

echo "# Done step2:      $(date)"

# ---------------------
# 03. Load mapped reads
# ---------------------

bam_file="$OUTDIR/bam/QUERY.adapter_filtered2.hg19_female.clean_srt_markdup.bam"
out_dir=$OUTDIR/Aneu_analysis
name="QUERY"
blacklist="sample_data/blacklist/hg19-blacklist-v1_id.bed"
genome_file="sample_data/genome_file/UCSC_hg19_female.fa.fai"

bash scripts/Step3_load_mapped_reads.sh \
$bam_file \
$out_dir \
$name $blacklist \
$genome_file \
&> "$OUTDIR/logs/Step3.log"

bash util/Fragment2bin.sh \
sample_data/fragment_data \
$OUTDIR/Aneu_analysis/bins/ \
sample_data/genome_file/UCSC_hg19_female.fa.fai \
&> "$OUTDIR/logs/Fragment2bin.log"

echo "# Done step3:      $(date)"

# ---------------------
# 04. MAD score
# ---------------------

bash scripts/Step4_MAD_score_QC.sh \
$OUTDIR/Aneu_analysis/bins/ \
$OUTDIR/Aneu_analysis/MAD_score/ \
test_run &> "$OUTDIR/logs/Step4.log"

echo "# Done step4:      $(date)"

# ---------------------
# 05-1. Check karyotype
# ---------------------

IDs=("P285_09_1" "P285_10_1" "P293_29_1" "P293_30_1" "P301_30_1")

for ID in ${IDs[@]}
do
bash scripts/Step5_1_Check_G1_cells.sh $OUTDIR/Aneu_analysis/bins/${ID}_mapq10_blacklist_bin.Rdata \
$OUTDIR/Aneu_analysis/G1_control/ &> $OUTDIR/logs/Step5_1_${ID}.log
done

echo "# Done step5-1:    $(date)"

# ---------------------
# 05-2. Merge G1
# ---------------------

out_merged_fragment="$OUTDIR/Aneu_analysis/G1_control/Merged_control_G1.fragment.Rdata"

bash scripts/Step5_2_Merge_G1_cells.sh \
$out_merged_fragment \
sample_data/fragment_data/P285_09_1_mapq10_blacklist_fragment.Rdata \
sample_data/fragment_data/P293_29_1_mapq10_blacklist_fragment.Rdata \
sample_data/fragment_data/P301_30_1_mapq10_blacklist_fragment.Rdata \
&> "$OUTDIR/logs/Step5_2.log"

echo "# Done step5-2:    $(date)"

# ---------------------
# Tag Density
# ---------------------

fragment_dir="sample_data/fragment_data"
outdir_tag="$OUTDIR/Aneu_analysis/TagDensity"
blacklist="sample_data/blacklist/hg19-blacklist-v1_id.bed"
genome_file="sample_data/genome_file/UCSC_hg19_female.fa.fai"

mkdir -p $outdir_tag
Rscript --vanilla util/TagDensity-w200ks40k.R $fragment_dir $outdir_tag $blacklist $genome_file &> "$OUTDIR/logs/TagDensity.log"

echo "# Done tagDensity: $(date)"

# ---------------------
# 06. Log2 med score
# ---------------------

fragment_file="sample_data/fragment_data/P293_29_1_mapq10_blacklist_fragment.Rdata"
out_dir="$OUTDIR/Aneu_analysis/Log2_Med"
ref_Rdata="$OUTDIR/Aneu_analysis/G1_control/Merged_control_G1.fragment.Rdata"
blacklist="sample_data/blacklist/hg19-blacklist-v1_id.bed"
genome_file="sample_data/genome_file/UCSC_hg19_female.fa.fai"

bash scripts/Step6_log2_median_RT_scores.sh \
$fragment_file $out_dir $ref_Rdata $blacklist $genome_file \
&> "$OUTDIR/logs/Step6.log"

echo "# Done step6:      $(date)"

# ---------------------
# 07. Log2 med score
# ---------------------

# 1-somy
binfile="$OUTDIR/Aneu_analysis/bins/P285_19_1_mapq10_blacklist_bin.Rdata"
somy="1-somy"
out_dir="$OUTDIR/Aneu_analysis/HMM/$somy"
ref_Rdata="$OUTDIR/Aneu_analysis/G1_control/Merged_control_G1.fragment.Rdata"
genome_file="sample_data/genome_file/UCSC_hg19_female.fa.fai"
binsize="80000"

bash scripts/Step7_Binarization.sh \
$binfile $out_dir $ref_Rdata $genome_file $binsize $somy \
&> "$OUTDIR/logs/Step7_1-somy.log"

# 2-somy
somy="2-somy"
out_dir="$OUTDIR/Aneu_analysis/HMM/$somy"

bash scripts/Step7_Binarization.sh \
$binfile $out_dir $ref_Rdata $genome_file $binsize $somy \
&> "$OUTDIR/logs/Step7_2-somy.log"

echo "# Done step7:      $(date)"

echo "# End   test run:  $(date)"
echo "# Test run finished successfully."

end=$(date +%s)
elapsed=$((end - start))

printf "##\n## Elapsed time: %02d min %02d sec\n##\n" $((elapsed / 60)) $((elapsed % 60))
