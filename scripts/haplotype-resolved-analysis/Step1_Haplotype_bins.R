#!/usr/bin/env Rscript

#v190401 Haplotype analysis Step1
#bin the bed files
#
#$ Rscript Step1_Haplotype_bins.R <blacklist> <genome_file> \
#  <allele_m_file> <allele_p_file> \
#  <allele_m_name> <allele_p_name> \
#  <outname> <bin_sizes(comma separated)>

library(AneuFinder)
options(scipen=100)
args=commandArgs(TRUE)

##loading black list and genome Info##
blacklist=args[1] #mm9-blacklist-id.bed
genome_tmp <- read.table(args[2],sep="\t") #UCSC_mm9.woYwR.fa.fai
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

allele_m_file=args[3] #cba bed file after rmdup
allele_p_file=args[4] #msm bed file after rmdup
allele_m_name=args[5] #cba
allele_p_name=args[6] #msm
outname=args[7] #Test

in_binsizes=args[8] #binsizes comma separated 400000,1000000
bin_sizes=strsplit(in_binsizes , ",")

haplo_bins_reads=list()
bins_allele_m <- binReads(paste(allele_m_file,sep=""), assembly=genome,
                     binsizes=as.numeric(bin_sizes[[1]]), chromosomes=chromosomes)
bins_allele_p <- binReads(paste(allele_p_file,sep=""), assembly=genome,
                     binsizes=as.numeric(bin_sizes[[1]]), chromosomes=chromosomes)

haplo_bins_reads[[allele_m_name]]=bins_allele_m
haplo_bins_reads[[allele_p_name]]=bins_allele_p
save(haplo_bins_reads,file=paste(outname,"_haplotype_bins.Rdata",sep=""))
