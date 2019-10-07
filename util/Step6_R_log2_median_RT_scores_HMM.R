#!/usr/bin/env Rscript
path="/usr/local/bin"
#scRepli-seq pipeline
#Computation of log2 med score from 2HMM bin file

library("AneuFinder")
source(normalizePath(paste0(path,"/util/Aneufinder_Optional_script.R")))
options(scipen=100)

args = commandArgs(TRUE)
HMMfile=args[1]
outdir=args[2]
#blacklist=args[3]
#genome_file=args[4]
#binsize=args[5]
#ref_Rdata=args[6]

median_scale = function (datalist) {
x=datalist
c=elementMetadata(x)$counts
med=median(c[c!=0])
c_m=log2(c/med)
out=data.frame(chr=seqnames(x),
                 starts=start(x)-1,
                 ends=end(x),
                 map_median_scale_log2=c_m)

return(out)
}


##loading black list and genome Info##
#genome_tmp <- read.table(genome_file,sep="\t") #
#genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
#chromosomes=as.character(genome$UCSC_seqlevel)

##load bin file
load(HMMfile)
test_map_med_log2=median_scale(binned.MAP.HMM.2$bins)
out=test_map_med_log2
a=out[order(out$chr,out$starts),]
a1=data.frame(chr=a[,1],s=a[,2],e=a[,3],c=a$map_median_scale_log2)
a2=a1[!is.infinite(a1$c),]

name=binned.MAP.HMM.2$ID
bsize=(a1[1,3] - a1[1,2])/1000

write.table(a2,paste(outdir,"/",name,"_",bsize,"k_2HMM_map_count_median_log2.bedGraph",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
