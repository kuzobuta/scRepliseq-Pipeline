#!/usr/bin/env Rscript
#Merge HMM Rdata into one txt file

library(AneuFinder)
options(scipen = 100)

args = commandArgs(TRUE)
HMMfile_dir=args[1]
outdir=args[2]
out_name=args[3]

HMMfiles=list.files(HMMfile_dir,pattern=".Rdata",full.names=T)

out=NULL;i=0
for (file in HMMfiles){
 i=i+1
 load(file) #HMM file
 name=binned.MAP.HMM.2$ID
 out_tmp=data.frame(chr=seqnames(binned.MAP.HMM.2[[2]]),
                    starts=start(binned.MAP.HMM.2[[2]])-1,
                    ends=end(binned.MAP.HMM.2[[2]]),
                    copy.numbers=elementMetadata(binned.MAP.HMM.2[[2]])$copy.number)
 if(i==1){
   out=out_tmp
   colnames(out)=c("chr","start","end",name)
 }else{
   out=cbind(out,out_tmp[,4])
   colnames(out)[dim(out)[2]]=name
 }
}

merge_out=out
write.table(merge_out,paste0(outdir,"/",out_name,"_merged_binary_data.txt"),sep="\t",col.names=T,row.names=F,quote=F)
