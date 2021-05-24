#!/usr/bin/env Rscript
#BrdUIP w200k s40k analysis using AneuFinder#
#From bam file
#Hisashi Miura <hisashi.miura[at]riken.jp>
#version: 210304
#
#$ Rscript <XXX.R> \
#<Early-S bam file> \
#<Late-S bam file> \
#<sample name> \
#<out directory> \
#<black list> \
#<genome file>
#
#This script generates .bedGraph file of BrdU-IP data 
#& one directories (fragment) under <out directory>
#Workflow:
#0. Load data from bam file [two files are saved as XXX_[EarlyS or LateS]_mapq10_blacklist_fragment.Rdata]
#1. Generate sliding window bins
#2. Load reads & rpm normalization
#3. E/(E+L) score scaled from -1 to 1
#4. Filtering the E/(E+L) score by <= 5% percentile of (E rpm + L rpm)

library("AneuFinder")
library("zoo")
library("pracma")

args = commandArgs(TRUE)
#fragment_E=args[1] #Early-S fraction
#fragment_L=args[2] #Late-S fraction
bam_E=args[1] #Early-S fraction (bam file)
bam_L=args[2] #Late-S fraction (bam file)
name=args[3]
outdir=args[4] 
blacklist=args[5] #e.x.) mm9.blacklist-id.bed
genome_file=args[6] #e.x.) 

dir.create(outdir,showWarnings=F)

options(scipen=100)

##loading black list and genome Info##
genome_tmp <- read.table(genome_file,sep="\t")
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)
##############################################################
##Extension of file name##
ext="_mapq10_blacklist_fragment.Rdata"
#ext2="_mapq10_blacklist_bin.Rdata"

##Setup output directories##
outdir_f=paste0(outdir,"/fragment")
dir.create(outdir,showWarnings = FALSE)
dir.create(outdir_f,showWarnings = FALSE)

#0. load data from bam file
##Load the bam file & save the fragment file (>10 MAPQ), filtering out the blacklist regions##
raw_reads=bam2GRanges(bam_E,remove.duplicate.reads = TRUE,min.mapq = 10,blacklist = blacklist)
save(raw_reads,file = paste0(outdir_f,"/",name,"_EarlyS_",ext))

raw_reads=bam2GRanges(bam_L,remove.duplicate.reads = TRUE,min.mapq = 10,blacklist = blacklist)
save(raw_reads,file = paste0(outdir_f,"/",name,"_LateS_",ext))

##############################################################
#1. Generate sliding window bins
#If you want to change the setting, please change these.
window=200000
sliding=40000
sw_name=paste0("w",round(window/1000),"ks",round(sliding/1000),"k")

chr_list=genome$UCSC_seqlevel
lcm=Lcm(window,sliding)
options(scipen=100)
print(paste("sliding size is",lcm/sliding,"for window:",window,"interval:",sliding))
sliding_window_data=list()
for (j in 1:(lcm/sliding)){
  tmp_out=NULL
  for (chr in chr_list){
    chr_size=genome$UCSC_seqlength[genome$UCSC_seqlevel==chr]
    s=(j-1)*sliding+window
    if (s>=chr_size){print(paste("skip",chr));next}
    tmp=data.frame(chr=chr,
                   start=seq(s,chr_size,window) - window + 1,
                   end=seq(s,chr_size,window))
    tmp_out=rbind(tmp_out,tmp)
  }
  sliding_window_data[[j]]=tmp_out
}

x_gr=list()
for (i in 1:length(sliding_window_data)){
  x_gr[[i]]=makeGRangesFromDataFrame(sliding_window_data[[i]])
}
##############################################################

##############################################################
norm_rpm_bins=function(x,rpm){
  out=NULL
  for (i in 1:(length(x)-1)){
    tmp=data.frame(chr=seqnames(x[[i]]),
                   starts=start(x[[i]])-1,
                   ends=end(x[[i]]),
                   rpm_counts=x[[i]]$counts * rpm)
    out=rbind(out,tmp)
  }
  out=out[order(out$starts),]
  srt=order(out$chr,out$starts)
  out2=out[srt,]
  out3=data.frame(chr=out2[,1],
                  starts=out2[,2] + (window - sliding)/2,
                  ends=out2[,2] + (window - sliding)/2 + sliding,
                  rpm_counts=out2[,4])
  return(out3)
}

norm_rpm_bins_fix=function(x,rpm){
  out=NULL
  tmp=data.frame(chr=seqnames(x[[1]]),
                 starts=start(x[[1]])-1,
                 ends=end(x[[1]]),
                 rpm_counts=x[[1]]$counts * rpm)
  srt=order(tmp$chr,tmp$starts)
  out2=tmp[srt,]
  out3=data.frame(chr=out2[,1],
                  starts=out2[,2],
                  ends=out2[,3],
                  rpm_counts=out2[,4])
  return(out3)
}

##############################################################
Percent_EL_IP=function(x){
  if(x[1]==0 && x[2]==0){
    score=NA
  }else if(x[1]==0 && x[2]>0){
    score=-1
  }else if (x[1]>0 && x[2]==0){
    score=1
  }else{
    score=(x[1]/(x[1]+x[2])-0.5)*2
  }
  return(score)
}
##############################################################

##############################################################
RT_IP_rpm_q_analysis = function(name,Efile,Lfile,genome,blacklist,x_gr){
  #2. Load reads & rpm normalization
  load(Efile)
  rpm=1000000/length(raw_reads)
  E_bins=binReads(raw_reads,
                  assembly=genome,
                  blacklist=blacklist,
                  bins=x_gr,
                  binsizes = 500000)
  E_bins_rpm=norm_rpm_bins(E_bins,rpm)
  load(Lfile)
  rpm=1000000/length(raw_reads)
  L_bins=binReads(raw_reads,
                  assembly=genome,
                  blacklist=blacklist,
                  bins=x_gr,
                  binsizes = 500000)
  L_bins_rpm=norm_rpm_bins(L_bins,rpm)
  #identical(E_bins_rpm[,1:3],L_bins_rpm[,1:3])
  
  x=cbind(E_bins_rpm,L_bins_rpm[,4])
  colnames(x)[4:5]=c("E_rpm","L_rpm")
  #3. E/(E+L) score scaled from -1 to 1
  Percent_E=apply(x[,4:5],1,Percent_EL_IP)
  IP_rpm_sw=x
  IP_rpm_sw$Percent_E=Percent_E
  IP_rpm_sw_c=IP_rpm_sw[complete.cases(IP_rpm_sw),]
  
  ##calcurate the percentile of "E + L rpm reads" w/o NA
  #q=quantile(IP_rpm_sw[,4]+IP_rpm_sw[,5],seq(0,1,0.01))
  EL_rpm = IP_rpm_sw_c[,4]+IP_rpm_sw_c[,5]
  q_c=quantile(EL_rpm[EL_rpm!=0],seq(0,1,0.01))
  
  #E_L_rpm=IP_rpm_sw[,4] + IP_rpm_sw[,5]
  E_L_rpm_c=IP_rpm_sw_c[,4] + IP_rpm_sw_c[,5] #E_rpm + L_rpm
  IP_rpm_sw_c_p=IP_rpm_sw_c
  
  ##Filter out 5 percentile data (E_rpm + L_rpm) or less 
  IP_rpm_sw_c_p[E_L_rpm_c<=q_c[6],6]=NA
  
  out=IP_rpm_sw_c_p[complete.cases(IP_rpm_sw_c_p),]
  write.table(out[,c(1,2,3,6)],paste0(outdir,"/",name,"_IP_Aneu_rmdup_",sw_name,"_Percent_q0.05.bedGraph"),
              sep="\t",col.names=F,row.names=F,quote=F)
}
##############################################################

##################Run###########################
fragment_E=paste0(outdir_f,"/",name,"_EarlyS_",ext)
fragment_L=paste0(outdir_f,"/",name,"_LateS_",ext)

RT_IP_rpm_q_analysis(name,fragment_E,fragment_L,genome,blacklist,x_gr)
