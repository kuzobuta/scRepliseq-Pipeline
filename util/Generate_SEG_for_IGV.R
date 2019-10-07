#!/usr/bin/env Rscript
path="/usr/local/bin"

#############################################################
#Generation of SEG file for visualizing multiple binarized scRepli-seq data
#In IGV, we showed the binarized scRepli-seq data in heatmap mode by the following scales; set Maximum scale to 1.0 for blue (replicated bins), set Minimum scale to -1.0 for yellow (unreplicated bins).
#Output file name should be name as "XXXX.seg"
#If you have a lot of samples (>100), the file size of "XXXX.seg" will be large (> 200Mb).
#

options(scipen = 100)
args = commandArgs(TRUE)
in_dirs=args[1] #bedGraph file directory
out_dir=args[2]
out_name=args[3]
sort_list_file=args[4]

################Get sorted file names#############
#sort_names need to be unique names
if (is.na(sort_list_file)){
  sort_names=NULL
}else{
  tmp=read.table(sort_list_file,sep="\t",header=T)
  sort_names=tmp[,1]
}

################Get file names####################
files=list.files(in_dirs,pattern = ".bedGraph")
files1=NULL
for (file in files){
  files1=c(files1,paste0(in_dirs,"/",file))
}

################Order the file names by the sorted order################
if (is.na(sort_list_file)){
  files1_srt=files1
}else{
  files1_srt=NULL
  for (n in sort_names){
    files1_srt=c(files1_srt,files1[grep(n,files1)])
  }
}

################load the bedGraph files and export into .seg file################

w=nchar(as.character(length(files1_srt)))
out=NULL;i=0
for (file in files1_srt){
  i=i+1
  out_i=formatC(i,width=w,flag="0")
  filename=paste0(out_i,"_",file)
  in_file=paste0(in_dirs,"/",file)
  if(i==1){
    tmp=read.table(in_file,sep="\t")
    tmp_seg.mean=tmp[,4]
    tmp_seg.mean[tmp[,4]==0]=NA
    out=cbind(filename,as.character(tmp[,1]),tmp[,2]+1,tmp[,3],tmp_seg.mean)
  }else{
    tmp=read.table(in_file,sep="\t")
    tmp_seg.mean=tmp[,4]
    tmp_seg.mean[tmp[,4]==0]=NA
    tmp_out=cbind(filename,as.character(tmp[,1]),tmp[,2]+1,tmp[,3],tmp_seg.mean)
    out=rbind(out,tmp_out)
  }
}
colnames(out)=c("'ID","chrom","loc.start","loc.end","seg.mean")
write.table(out,paste0(out_dir,"/",out_name,".seg"),
            sep="\t",col.names=T,row.names=F,quote=F)

