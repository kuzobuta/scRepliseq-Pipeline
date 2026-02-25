#!/usr/bin/env Rscript
#fragment file to bin files by AneuFinder
#Hisashi Miura <hisashi.miura@riken.jp>
#
#$ Rscript <XXX.R> <fragment file directory> <out_dir> <genome file> <bin size> 
#<bin size> is a option
#If <bin size> is not specified, default binsizes are selected as 40,80,100,200,500kb & 1Mb
#Input fragment R file is XXX_mapq10_blacklist_fragment.Rdata

args = commandArgs(TRUE)
frag_dir=args[1]
out_dir=args[2]
genome_file=args[3]
binsize=args[4]
options(scipen=100)

if(!is.na(binsize)){
 binsizes=binsize
}else{
 binsizes=c(40000,80000,100000,200000,500000,1000000)
}

##Extension of file name##
ext="_mapq10_blacklist_fragment.Rdata"
ext2="_mapq10_blacklist_bin.Rdata"

library(AneuFinder)

##loading black list and genome Info##
genome_tmp <- read.table(genome_file,sep="\t") #UCSC_mm9.woYwR.fa.fai
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

files=list.files(frag_dir,pattern=ext,full.names=T)

for (fragfile in files){
 ##load the fragment R file including >10 MAPQ reads & reads excluding from the blacklist regions##
 load(fragfile) #as raw_reads

 name=basename(fragfile)
 name=sub(ext, "", name)
 
 ##save the bin data file ##
 bins_reads=binReads(raw_reads,
                     assembly=genome,
                     chromosomes=chromosomes,
                     binsizes=binsizes)
 names(bins_reads) = binsizes
 rpm=1000000/length(raw_reads)
 bins_reads[["rpm"]]=rpm
 save(bins_reads,file=paste(out_dir,"/",name,ext2,sep=""))
}
