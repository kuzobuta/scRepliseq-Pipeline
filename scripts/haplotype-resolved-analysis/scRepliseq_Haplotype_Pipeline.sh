#!/bin/bash
#v201124
#Example scripts for haplotype-resolved scRepli-seq analysis
#Please use "AneuFinder 1.2.1"

###################cat gz file##############
cat sample_data/CBMS1_d0_G1_Merged.cba.rmdup.bed.gz-* > sample_data/CBMS1_d0_G1_Merged.cba.rmdup.bed.gz 
cat sample_data/CBMS1_d0_G1_Merged.msm.rmdup.bed.gz-* > sample_data/CBMS1_d0_G1_Merged.msm.rmdup.bed.gz 

###################Step1 load maternal/paternal allele bed file##################
resultdir="test_results"
blacklist="sample_data/mm9-blacklist-id.bed"
genome_file="sample_data/UCSC_mm9.woYwR.fa.fai"
allele_m_name="cba"
allele_p_name="msm"

mkdir -p ${resultdir}

###################For control G1################
#d0
allele_m_file_G1="sample_data/CBMS1_d0_G1_Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_G1="sample_data/CBMS1_d0_G1_Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_G1="d0_G1"
bin_sizes="400000,1000000"
outdir=${resultdir}/step1_Control_G1

mkdir -p ${outdir}
Rscript Step1_Haplotype_bins.R ${blacklist} ${genome_file} \
${allele_m_file_G1} ${allele_p_file_G1} \
${allele_m_name} ${allele_p_name} \
${outdir}/${outname_G1} \
${bin_sizes}

#d7
allele_m_file_G1="sample_data/CBMS1_d7_G1_Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_G1="sample_data/CBMS1_d7_G1_Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_G1="d7_G1"
bin_sizes="400000,1000000"
outdir=${resultdir}/step1_Control_G1

Rscript Step1_Haplotype_bins.R ${blacklist} ${genome_file} \
${allele_m_file_G1} ${allele_p_file_G1} \
${allele_m_name} ${allele_p_name} \
${outdir}/${outname_G1} \
${bin_sizes}

#Outputfile (1file):
#${outname}_haplotype_bins.Rdata #as Rdata, name: haplo_bins_reads

###################For S-phase################
###For only 1 file###
bin_sizes="400000,1000000"
outdir=${resultdir}/step1_S_phase
mkdir -p ${outdir}

#d0
allele_m_file_S="sample_data/d0_MidS_P263_11_1.Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_S="sample_data/d0_MidS_P263_11_1.Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_S="d0_MidS_P263_11"

Rscript Step1_Haplotype_bins.R ${blacklist} ${genome_file} \
${allele_m_file_S} ${allele_p_file_S} \
${allele_m_name} ${allele_p_name} \
${outdir}/${outname_S} \
${bin_sizes}

#d7
allele_m_file_S="sample_data/d7_MidS_P293_79_1_R1_A_SEQXE_v2.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_S="sample_data/d7_MidS_P293_79_1_R1_A_SEQXE_v2.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_S="d7_MidS_P293_79"

Rscript Step1_Haplotype_bins.R ${blacklist} ${genome_file} \
${allele_m_file_S} ${allele_p_file_S} \
${allele_m_name} ${allele_p_name} \
${outdir}/${outname_S} \
${bin_sizes}


#Outputfile (1file):
#${outname}_haplotype_bins.Rdata #as Rdata, name: haplo_bins_reads

###For loop S-phase samples###
#
#in_dir=""
#bin_sizes="400000,1000000"
#outdir=${resultdir}/step1_S_phase
#mkdir -p ${outdir}

#in_dir="" #input file directory (allele maternal & paternal bed files)
#files=(`ls ${in_dir}/*.${allele_m_name}.rmdup.bed.gz`)
#for file in ${files[@]}
#do
#tmp_mat_file=${file} #XXX.cba.rmdup.bed.gz$
#prefix=`basename $file`
#prefix=${prefix%.${allele_m_name}.rmdup.bed.gz}
#tmp_pat_file=${in_dir}/${prefix}.${allele_p_name}.rmdup.bed.gz #XXX.msm.rmdup.bed.gz$
#Rscript Step1_Haplotype_bins.R ${blacklist} ${genome_file} \
#${tmp_mat_file} ${tmp_pat_file} \
#${allele_m_name} ${allele_p_name} \
#${outdir}/${prefix} \
#${bin_sizes}
#done

#Outputfile (1file):
#${prefix}_haplotype_bins.Rdata #as Rdata, name: haplo_bins_reads
#in ${outdir}

###################Step2.1 Check copy numbers in maternal/paternal allele merged G1 bed file##################
#To generate the merged G1 bed file
#$ cat XXX.cba.bed.gz YYY.cba.bed.gz > ZZZ.cba.merged.bed.gz
#HMM analysis was fixed in 400kb

outdir=${resultdir}/step2_karyogram
mkdir -p ${outdir}
bin_sizes="400000,1000000"

allele_m_file_G1="sample_data/CBMS1_d0_G1_Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_G1="sample_data/CBMS1_d0_G1_Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_G1="d0_G1"

Rscript Step2_Haplotype_Merged_G1_HMM_check.R ${blacklist} ${genome_file} \
${allele_m_file_G1} ${allele_p_file_G1} \
${allele_m_name} ${allele_p_name} \
${outdir}/${outname_G1}

allele_m_file_G1="sample_data/CBMS1_d7_G1_Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_G1="sample_data/CBMS1_d7_G1_Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_G1="d7_G1"

Rscript Step2_Haplotype_Merged_G1_HMM_check.R ${blacklist} ${genome_file} \
${allele_m_file_G1} ${allele_p_file_G1} \
${allele_m_name} ${allele_p_name} \
${outdir}/${outname_G1}

#Outputfiles (3files):
#${outname_G1}_karyogram_merged_G1_${allele_m_name}.pdf #karyogram plot by copynumber analysis
#${outname_G1}_karyogram_merged_G1_${allele_p_name}.pdf #karyogram plot by copynumber analysis
#${outname_G1}_${allele_m}_${allele_p}_G1_control_400k_4HMM.Rdata #as Rdata, name: bins_400k_HMM_haplo

###################Step2.2 Select the 2 copy chr or regions##################
#Based on the results of step2.1, please select the 2 copy chr or positions and filter out regions
#This result may be differ among the cell batches

###################Step3 Generate custom blacklist for each allele##################
# The genomic bins with > Med+1.5*IQR or < Med-1.5*IQR are defined as a custom blacklist,
# & filtered out for further analysis. 
# ChrX is separately computed as the read coverage is differ from autosome.
# As an input file, you need to generate the "discriminated_regions.bed" file.
# If you want to select all chr, you can set NA in <start> / <end>.
#
# "high_copy"  : high copy region [Med +/- 1.5*IQR are separately computed] 
# "filter_out" : filtering out regions
# "X_specific" : X chromosome specific options
#
# Example of "discriminated_regions.bed" file for d7
# <chr_name>  <start>	<end>       <allele_name>   <type>
# chr8	      NA        NA          cba             high_copy
# chr12       74800000  121200000   cba             high_copy
# chr18       NA        NA          cba             filter_out
# chr18       NA        NA          msm             filter_out
# chrX        NA        NA          cba             X_specific
# chrX        NA        NA          msm             X_specific
#
# For filtering out chrX for further analysis [this might be used for male samples]
# chrX        NA		NA			cba				filter_out
# chrX        NA		NA			msm				filter_out
#
# If all chromosomes are normal and samples are female, you may use the following example
# chrX		  NA		NA			cba				X_specific
# chrX		  NA		NA			msm				X_specific
# 

bin_size=400000
outdir=${resultdir}/step3_blacklist
mkdir -p ${outdir}

#d0_G1
allele_m_file_G1="sample_data/CBMS1_d0_G1_Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_G1="sample_data/CBMS1_d0_G1_Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_G1="d0_G1"
dis_bed_file="sample_data/Test_d0_G1_discriminated_regions.bed"

Rscript Step3_Haplotype_blacklist.R ${blacklist} ${genome_file} \
${allele_m_file_G1} ${allele_p_file_G1} \
${allele_m_name} ${allele_p_name} \
${bin_size} \
${dis_bed_file} \
${outdir}/${outname_G1}

#d7_G1
allele_m_file_G1="sample_data/CBMS1_d7_G1_Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_G1="sample_data/CBMS1_d7_G1_Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_G1="d7_G1"
dis_bed_file="sample_data/Test_d7_G1_discriminated_regions.bed"

Rscript Step3_Haplotype_blacklist.R ${blacklist} ${genome_file} \
${allele_m_file_G1} ${allele_p_file_G1} \
${allele_m_name} ${allele_p_name} \
${bin_size} \
${dis_bed_file} \
${outdir}/${outname_G1}

#Outputfiles:
#${outname_G1}_XX_${allele_name}_Type_IQRfilter.pdf #Plot for blacklist 
#${outname_G1}_${allele_m_name}_${allele_p_name}_blacklist.Rdata #R data as "blacklist_haplotype"
#${outname_G1}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata #R data name as "haplo_control_S_G1_reads_newBL"

###################Step4 allele specific copy number analysis using HMM##################
bin_size=400000
outdir=${resultdir}/step4_HMM
mkdir -p ${outdir}

##d0 MidS sample##
R_haplotype_file="${resultdir}/step1_S_phase/d0_MidS_P263_11_haplotype_bins.Rdata"
outname_G1="d0_G1"
somy="2-somy" #1-somy or 2-somy

Ref_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata"
black_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_blacklist.Rdata"
Rscript Step4_Haplotype_HMM_somy_mode_different_binsizes.R ${genome_file} \
${R_haplotype_file} \
${allele_m_name} ${allele_p_name} \
${bin_size} \
${outdir} \
${Ref_file} \
${black_file} \
${somy}

##d7 MidS sample##
R_haplotype_file="${resultdir}/step1_S_phase/d7_MidS_P293_79_haplotype_bins.Rdata"
outname_G1="d7_G1"
somy="2-somy" #1-somy or 2-somy

Ref_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata"
black_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_blacklist.Rdata"
Rscript Step4_Haplotype_HMM_somy_mode_different_binsizes.R ${genome_file} \
${R_haplotype_file} \
${allele_m_name} ${allele_p_name} \
${bin_size} \
${outdir} \
${Ref_file} \
${black_file} \
${somy}

#Outputfiles:
#It generates 
#"Rdata"  : HMM results scored in Rdata
#"binary" : bedGraph format
#"plot"   : 2HMM read distribution plot
#directories
#

###For loop S-phase samples###
#
#bin_size=400000
#outdir=${resultdir}/step4_HMM
#somy="2-somy" #1-somy or 2-somy
#mkdir -p ${outdir}
#
#outname_G1=""
#
#Ref_file=${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata
#black_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_blacklist.Rdata"
#
#in_dir=${resultdir}/step1_S_phase
#files(`ls ${in_dir}/*_haplotype_bins.Rdata`)
#
#for file in ${files[@]}
#do
#f=`basename $file`
#Sname=${f%_haplotype_bins.Rdata}
#R_haplotype_file=${resultdir}/step1_S_phase/${Sname}_haplotype_bins.Rdata
#Rscript Step4_Haplotype_HMM_somy_mode_different_binsizes.R ${genome_file} \
#${R_haplotype_file} \
#${allele_m_name} ${allele_p_name} \
#${bin_size} \
#${outdir} \
#${Ref_file} \
#${black_file} \
#${somy}
#done
#

###################Step5 allele specific log2 median (RT score) by w1Ms40k ##################
#This analysis will take ~5min per sample..
outdir=${resultdir}/step5_log2med
mkdir -p ${outdir}

##d0 MidS sample##
allele_m_file_S="sample_data/d0_MidS_P263_11_1.Merged.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_S="sample_data/d0_MidS_P263_11_1.Merged.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_S="d0_MidS_P263_11"
outname_G1="d0_G1"

Ref_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata"
black_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_blacklist.Rdata"
Rscript Step5_Haplotype_log2med.R ${blacklist} ${genome_file} \
${allele_m_file_S} ${allele_p_file_S} \
${allele_m_name} ${allele_p_name} \
${outdir} \
${Ref_file} \
${black_file} \
${outname_S}

##d7 MidS sample##
allele_m_file_S="sample_data/d7_MidS_P293_79_1_R1_A_SEQXE_v2.cba.rmdup.bed.gz" #XXX.cba.rmdup.bed.gz$
allele_p_file_S="sample_data/d7_MidS_P293_79_1_R1_A_SEQXE_v2.msm.rmdup.bed.gz" #XXX.msm.rmdup.bed.gz$
outname_S="d7_MidS_P293_79"
outname_G1="d7_G1"

Ref_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata"
black_file="${resultdir}/step3_blacklist/${outname_G1}_${allele_m_name}_${allele_p_name}_blacklist.Rdata"
Rscript Step5_Haplotype_log2med.R ${blacklist} ${genome_file} \
${allele_m_file_S} ${allele_p_file_S} \
${allele_m_name} ${allele_p_name} \
${outdir} \
${Ref_file} \
${black_file} \
${outname_S}



