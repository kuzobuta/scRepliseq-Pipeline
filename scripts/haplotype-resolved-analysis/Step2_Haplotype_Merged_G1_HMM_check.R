#!/usr/bin/env Rscript

#v190401 Haplotype analysis Step2
#HMM using merged G1 cells
#
#$ Rscript Step2_Haplotype_Merged_G1_HMM_check.R \
#  <blacklist> <genome_file> \
#  <Merged_G1_allele_m> <Merged_G1_allele_p> \
#  <allele_m_name> <outname>

library(AneuFinder)
library(GenomicRanges)
library(ggplot2)
options(scipen=100)
args=commandArgs(TRUE)

##loading black list and genome Info##
blacklist=args[1] #mm9-blacklist-id.bed
genome_tmp <- read.table(args[2],sep="\t") #UCSC_mm9.woYwR.fa.fai
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

allele_m_bedfile=args[3] #cba bed file
allele_p_bedfile=args[4] #msm bed file

allele_m_name=args[5] #cba
allele_p_name=args[6] #msm

outname=args[7] #


#######################################cba & msm HMM#######################################

bins_400k_HMM_haplo=list()
bins_allele_m <- binReads(allele_m_bedfile, assembly=genome,
                 binsizes=400000, chromosomes=chromosomes)

bins_400k_HMM_haplo[[allele_m_name]]=findCNVs(bins_allele_m$`400000`,
                                      method = 'HMM',num.threads=10,
                                      eps=0.01,
                                      states=c("zero-inflation","0-somy","1-somy","2-somy","3-somy"),
                                      max.iter = 3000)

pdf(paste0(outname,"_karyogram_merged_G1_",allele_m_name,".pdf"))
plot(bins_400k_HMM_haplo$cba,type="karyogram")
dev.off()

bins_allele_p <- binReads(allele_p_bedfile, assembly=genome,
                 binsizes=400000, chromosomes=chromosomes)
bins_400k_HMM_haplo[[allele_p_name]]=findCNVs(bins_allele_p$`400000`,
                                      method = 'HMM',num.threads=10,
                                      eps=0.01,
                                      states=c("zero-inflation","0-somy","1-somy","2-somy","3-somy"),
                                      max.iter = 3000)

pdf(paste0(outname,"_karyogram_merged_G1_",allele_p_name,".pdf"))
plot(bins_400k_HMM_haplo$msm,type="karyogram")
dev.off()

save(bins_400k_HMM_haplo,file=paste0(outname,"_",allele_m_name,"_",allele_p_name,"_G1_control_400k_4HMM.Rdata"))
