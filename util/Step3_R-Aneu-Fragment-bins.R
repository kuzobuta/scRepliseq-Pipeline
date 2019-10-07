#!/usr/bin/env Rscript
#scRepli-seq pipeline
#Importing the bam/bed file in Rdata by AneuFinder 
#fragment and bin analysis
#Hisashi Miura <hisashi.miura@riken.jp>
#
#$ Rscript <XXX.R> <bam/bed file> <out directory> <sample name> <black list> <genome file>
#This script generates two directories (fragment, bins) under <out directory>

args = commandArgs(TRUE)
bamfile=args[1]
out_dir=args[2]
name=args[3]
blacklist=args[4]
genome_file=args[5]
options(scipen=100)

##Extension of file name##
ext="_mapq10_blacklist_fragment.Rdata"
ext2="_mapq10_blacklist_bin.Rdata"

library(AneuFinder)

##loading black list and genome Info##
genome_tmp <- read.table(genome_file,sep="\t") #UCSC_mm9.woYwR.fa.fai
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

##setup output directories##
out_dir_f=paste0(out_dir,"/fragment")
out_dir_b=paste0(out_dir,"/bins")
dir.create(out_dir,showWarnings = FALSE)
dir.create(out_dir_f,showWarnings = FALSE)
dir.create(out_dir_b,showWarnings = FALSE)

##save the fragment file (>10 MAPQ), filtering out the blacklist regions##
raw_reads=bam2GRanges(bamfile,remove.duplicate.reads = TRUE,min.mapq = 10,blacklist = blacklist)
save(raw_reads,file = paste0(out_dir_f,"/",name,ext))

##save the bin data file ##
bins_reads=binReads(raw_reads,
                    assembly=genome,
                    chromosomes=chromosomes,
                    binsizes=c(40000,80000,100000,200000,500000))

rpm=1000000/length(raw_reads)
bins_reads[["rpm"]]=rpm
save(bins_reads,file=paste(out_dir_b,"/",name,ext2,sep=""))
