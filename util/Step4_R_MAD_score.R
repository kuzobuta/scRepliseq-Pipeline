#!/usr/bin/env Rscript
#scRepli-seq pipeline

library(AneuFinder)

options(scipen=100)
args = commandArgs(TRUE)
in_dir=args[1]
out_dir=args[2]
out_name=args[3]

files=list.files(in_dir,pattern=".Rdata",full.names=T)

out=NULL
for (file in files){
load(file) #bin file

#200k bins data
data=data.frame(chr=seqnames(bins_reads$`200000`),
                starts=start(bins_reads$`200000`)-1,
                ends=end(bins_reads$`200000`),
                elementMetadata(bins_reads$`200000`)$counts)

rmZero=data[,4] !=0
data_rmZero=data[rmZero,]
med=median(data_rmZero[,4])
med_scaled=data_rmZero[,4]/med
MAD=mad(med_scaled)
MAD_log=mad(log2(med_scaled))

name=basename(file)

out=rbind(out,c(name,MAD_log))
}

colnames(out)=c("Name","MAD_score_log2")

write.table(out,paste0(out_dir,"/",out_name,"_MAD_scores_log2.txt"),sep="\t",row.names=F,col.names=T,quote=F)
