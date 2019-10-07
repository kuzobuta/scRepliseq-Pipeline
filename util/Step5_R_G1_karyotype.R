#!/usr/bin/env Rscript
#6HMM by 500kb for checking the karyotype

library(AneuFinder)

args = commandArgs(TRUE)
binfile=args[1]
out_dir=args[2]

name=basename(binfile)
name=sub("_mapq10_blacklist_bin.Rdata","",name)

load(binfile)
bins_500k_HMM=findCNVs(bins_reads$`500000`,eps=0.01,
                              method = 'HMM',
                              states=c("zero-inflation","0-somy","1-somy","2-somy","3-somy","4-somy","5-somy","6-somy"),
                              max.iter = 3000)

bins_500k_HMM$ID=c(name)
save(bins_500k_HMM,file=paste(out_dir,"/",name,"_500k_6HMM.Rdata",sep=""))

name=bins_500k_HMM$ID
pdf(paste(out_dir,"/",name,"_500k_6HMM_karyogram.pdf",sep=""))
print(plot(bins_500k_HMM,type="karyogram"))
dev.off()
