#!/usr/bin/env Rscript
#Repliscores of the inpuf samples using binarized scRepli-seq data

library(AneuFinder)
options(scipen = 100)

args = commandArgs(TRUE)
HMMfile_dir=args[1]
out_dir=args[2]
out_name=args[3]
order_rep=args[4] #1 for order based on the repliscores, 0 is not.

HMMfiles=list.files(HMMfile_dir,pattern=".Rdata",full.names=T)

out=NULL
for (file in HMMfiles){
 load(file) #HMM file
 name=binned.MAP.HMM.2$ID
 tmp=binned.MAP.HMM.2[[2]]
 c=tmp$copy.number
 repscore_total=sum(c==2)/sum(c!=0)*100
 
 tmp_woX=tmp[seqnames(tmp)!="chrX"]
 c=tmp_woX$copy.number
 repscore_woX=sum(c==2)/sum(c!=0)*100
 
 out=rbind(out,c(name,repscore_total,repscore_woX))
}

colnames(out)=c("Name","Repliscore_Total","Repliscore_without_X")

if(order_rep==1){
 out_srt=out[order(out[,3]),]
}else{
 out_srt=out
}

write.table(out_srt,paste0(out_dir,"/",out_name,"_Repliscores.txt"),sep="\t",row.names=F,col.names=T,quote=F)
