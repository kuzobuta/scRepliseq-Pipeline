#!/usr/bin/env Rscript

#v190403 Haplotype log2med (RT score) analysis Step5
#
#
#$ Rscript Step5_Haplotype_log2med.R <blacklist> <genome_file> \
#  <allele_m_file> <allele_p_file> \
#  <allele_m_name> <allele_p_name> \
#  <outdir> \
#  <Ref_file> \
#  <black_file> \
#  <name>

library("zoo")
library("pracma")
library("AneuFinder")
source("Aneufinder_Optional_script.R")
args = commandArgs(TRUE)

#####loading black list and genome Info####
blacklist=args[1] #mm9-blacklist-id.bed
genome_tmp <- read.table(args[2],sep="\t") #UCSC_mm9.woYwR.fa.fai
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

allele_m_file=args[3] #cba bed file after rmdup
allele_p_file=args[4] #msm bed file after rmdup
allele_m_name=args[5] #cba
allele_p_name=args[6] #msm

out_dir=args[7] #Test
Ref_file=args[8] #From step 3 as haplo_control_S_G1_reads_newBL : ${outname}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata
black_file=args[9] #blacklist_haplotype : ${outname}_${allele_m_name}_${allele_p_name}_blacklist.Rdata

name=args[10]

out_dir_allele_m=paste(out_dir,"/",allele_m_name,sep="")
out_dir_allele_p=paste(out_dir,"/",allele_p_name,sep="")
dir.create(out_dir_allele_m,showWarnings = F)
dir.create(out_dir_allele_p,showWarnings = F)

load(Ref_file)   #blacklist filtered reads
load(black_file) #blacklist_haplotype

#############
window=1000000
#sliding=40000
sliding=200000 #change the condition as Nat Protoc. 
W="1M"
#S="40k"
S="200k" #change the condition as Nat Protoc.

####Sliding window function####
lcm=Lcm(window,sliding)
options(scipen=100)
print(paste("sliding size is",lcm/sliding,"for window:",window,"interval:",sliding))
sliding_window_data=list()
for (j in 1:(lcm/sliding)){
  tmp_out=NULL
  for (chr in chromosomes){
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

####Median function####
median_scale=function(datalist){
  x=datalist[[1]]
  c=elementMetadata(x)$counts
  med=median(c[c!=0])
  c_m=log2(c/med)
  out=data.frame(chr=seqnames(x),
                 starts=start(x)-1,
                 ends=end(x),
                 map_median_scale_log2=c_m)
  return(out)
}

####load data as reads####
#Bed file
allele_m_reads=bed2GRanges(allele_m_file,assembly = genome,chromosomes = chromosomes,blacklist = blacklist)
allele_p_reads=bed2GRanges(allele_p_file,assembly = genome,chromosomes = chromosomes,blacklist = blacklist)
#Bam file
#allele_m_reads=bam2GRanges(allele_m_file,assembly = genome,chromosomes = chromosomes,blacklist = blacklist)
#allele_p_reads=bam2GRanges(allele_p_file,assembly = genome,chromosomes = chromosomes,blacklist = blacklist)

allele_m_bins=binReads(allele_m_reads,
                       assembly=genome,
                       blacklist=blacklist_haplotype[[allele_m_name]],
                       bins=x_gr)

allele_p_bins=binReads(allele_p_reads,
                       assembly=genome,
                       blacklist=blacklist_haplotype[[allele_p_name]],
                       bins=x_gr)

####Mappability Correction####
options(warn=-1)
sliding_window_map=list()
for (i in 1:(length(allele_m_bins)-1)){
  ##Filter the blacklist reads (read count to zero)##
  bins=allele_m_bins[[i]]
  ov=findOverlaps(bins,blacklist_haplotype[[allele_m_name]])
  bins_filter=bins
  bins_filter$counts[queryHits(ov)]=0
  bins_filter$pcounts[queryHits(ov)]=0
  bins_filter$mcounts[queryHits(ov)]=0
  sliding_window_map[[allele_m_name]][[i]]=correctMappability(bins_filter,
                                                              same.binsize = T,
                                                              haplo_control_S_G1_reads_newBL[[allele_m_name]],
                                                              remove.duplicate.reads=F,
                                                              assembly=genome)
  bins=allele_p_bins[[i]]
  ov=findOverlaps(bins,blacklist_haplotype[[allele_p_name]])
  bins_filter=bins
  bins_filter$counts[queryHits(ov)]=0
  bins_filter$pcounts[queryHits(ov)]=0
  bins_filter$mcounts[queryHits(ov)]=0
  sliding_window_map[[allele_p_name]][[i]]=correctMappability(bins_filter,
                                                              same.binsize = T,
                                                              haplo_control_S_G1_reads_newBL[[allele_p_name]],
                                                              remove.duplicate.reads=F,
                                                              assembly=genome)
}

####Log2 median####
#allele_m
map_med_log2=NULL
for (i in 1:(length(allele_m_bins)-1) ){
  map_med_log2=rbind(map_med_log2,median_scale(sliding_window_map[[allele_m_name]][[i]]))
}
out=map_med_log2
a=out[order(out$chr,out$starts),]
a1=data.frame(chr=a[,1],s=round((a[,2]+a[,3])/2-sliding/2),e=round((a[,2]+a[,3])/2-sliding/2)+sliding,c=a$map_median_scale_log2)
a2=a1[!is.infinite(a1$c),]
write.table(a2,paste(out_dir_allele_m,"/",name,"_",allele_m_name,"_w",W,"s",S,"_map_count_median_log2.bedGraph",sep=""),col.names=F,row.names=F,sep="\t",quote=F)

#allele_p
map_med_log2=NULL
for (i in 1:(length(allele_m_bins)-1) ){
  map_med_log2=rbind(map_med_log2,median_scale(sliding_window_map[[allele_p_name]][[i]]))
}
out=map_med_log2
a=out[order(out$chr,out$starts),]
a1=data.frame(chr=a[,1],s=round((a[,2]+a[,3])/2-sliding/2),e=round((a[,2]+a[,3])/2-sliding/2)+sliding,c=a$map_median_scale_log2)
a2=a1[!is.infinite(a1$c),]
write.table(a2,paste(out_dir_allele_p,"/",name,"_",allele_p_name,"_w",W,"s",S,"_map_count_median_log2.bedGraph",sep=""),col.names=F,row.names=F,sep="\t",quote=F)



