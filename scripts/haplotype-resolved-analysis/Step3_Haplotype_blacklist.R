#!/usr/bin/env Rscript

#v190401 Haplotype analysis Step2.2
#Generation of the blacklist from HMM results
#
#$ Rscript Step3_Haplotype_blacklist.R <blacklist> <genome_file> \
#  <allele_m_bedfile> <allele_p_bedfile> \
#  <allele_m_name> <allele_p_name> \
#  <bin_size> <discriminated_bed_file> \
#  <outname>
#
#For calculating med +- 1.5 * IQR, chrX is computed separately from autosome,
#as the read coverage of chrX is slightly lower than the read coverage of autosome.
#Example of default bed file : "Test_default_discriminated_regions.bed"
#As an input file, you need to generate the "discriminated_regions.bed" file.
#If you want to select all chr, you can set NA in <start> / <end>.
#type :
#"high_copy"  : high copy number region, calculate the med +- 1.5 * IQR within this region [Optional]
#"filter_out" : filtering out region [Optional, you may use this option for filtering out abnormal chromosomes]
#"X_specific" : default
#Example of "discriminated_regions.bed" file
#<chr_name>	<start>	<end>		<allele_name>	<type>
#chr8		NA		NA			cba				high_copy
#chr12		7480000	121200000	cba				high_copy
#chr18		NA		NA			cba				filter_out
#chr18		NA		NA			msm				filter_out
#chrX		NA		NA			cba				X_specific
#chrX		NA		NA			msm				X_specific

library(AneuFinder)
library(GenomicRanges)
library(ggplot2)
options(scipen=100)
args=commandArgs(TRUE)

allele_m_bedfile=args[3]
allele_p_bedfile=args[4]

allele_m_name=args[5] #cba
allele_p_name=args[6] #msm

bin_size=args[7] #400000

dis_bedfile=args[8] #discriminated bed file

outname=args[9] #

########functions######
calc_black_list=function(bins_data,dis_regions,allele_name){
  blacklist=NULL
  x=bins_data
  tmp=dis_regions[dis_regions[,4]==allele_name,]
  #Set NA to entire chromosome
  for (i in 1:dim(tmp)[1]){
    tmp2=tmp[i,]
    if(is.na(tmp2[2])){
      start=1
      end=genome[which(genome$UCSC_seqlevel==as.character(tmp2[1,1])),2]
      tmp[i,2]=start - 1
      tmp[i,3]=end
    }
  }
  #Filterout genomic regions from bins_data
  fo_regions=tmp[tmp[,5]=="filter_out",]
  if(dim(fo_regions)[1]>0){
    fo_gr=GRanges(seqnames=fo_regions[,1],IRanges(start=fo_regions[,2]+1,end=fo_regions[,3]))
    y=subsetByOverlaps(x,fo_gr,invert=T)
    blacklist_tmp = granges(fo_gr)
    if (is.null(blacklist)){
      blacklist=blacklist_tmp
    }else{
      blacklist=c(blacklist,blacklist_tmp)
    }
  }else{
    y=x
  }
  
  #Compute "Med +- 1.5 * IQR" in the discriminated regions (e.x. high_copy, X_specific regions)
  j=0
  for (i in 1:dim(tmp)[1]){
    if(tmp[i,5]=="filter_out"){next}else{
      j=j+1
      #Select the specific genomic regions set in bed file
      gr=GRanges(seqnames=tmp[i,1],IRanges(start=tmp[i,2]+1,end=tmp[i,3]))
      z=subsetByOverlaps(y,gr)
      #Med +- 1.5xIQR
      c=z$counts
      c_c=c[c>0]
      iqr=IQR(c_c)
      l=median(c_c) - 1.5*iqr
      u=median(c_c) + 1.5*iqr
      #Plot the threshold lines as red (Med + 1.5 x IQR, Med - 1.5 x IQR)
      pdf(paste0(outname,"_",j,"_",allele_name,"_",as.character(tmp[i,1]),"_",as.character(tmp[i,5]),"_IQRfilter.pdf")) #filtering plot
      p <- plot(z)
      p <- p + geom_hline(aes(yintercept=l), color='red')
      p <- p + geom_hline(aes(yintercept=u), color='red')
      print(p)
      dev.off()
      rm(p)
      blacklist_tmp=granges(z[z$counts > u | z$counts < l])
      if (is.null(blacklist)){
        blacklist=blacklist_tmp
      }else{
        blacklist=c(blacklist,blacklist_tmp)}
    }
  }
  
  #Compute "Med +- 1.5 * IQR" in the rest of chromosomes
  gr=GRanges(seqnames=tmp[,1],IRanges(start=tmp[,2]+1,end=tmp[,3]))
  z=subsetByOverlaps(y,gr,invert=T)
  c=z$counts
  c_c=c[c>0]
  iqr=IQR(c_c)
  l=median(c_c) - 1.5*iqr
  u=median(c_c) + 1.5*iqr
  pdf(paste0(outname,"_",j+1,"_",allele_name,"_normal_copies_IQRfilter.pdf")) #filtering plot
  p <- plot(z)
  p <- p + geom_hline(aes(yintercept=l), color='red')
  p <- p + geom_hline(aes(yintercept=u), color='red')
  print(p)
  dev.off()
  blacklist_tmp=granges(z[z$counts > u | z$counts < l])
  if (is.null(blacklist)){
    blacklist=blacklist_tmp
  }else{
    blacklist=c(blacklist,blacklist_tmp)}
  return(blacklist)
}

###loading black list and genome Info###
blacklist=args[1]
genome_tmp <- read.table(args[2],sep="\t") #UCSC_mm9.woYwR.fa.fai
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

#load allele maternal
bins_allele_m <- binReads(allele_m_bedfile, assembly=genome,binsizes=as.numeric(bin_size),chromosomes = chromosomes)

#load allele paternal
bins_allele_p <- binReads(allele_p_bedfile, assembly=genome,binsizes=as.numeric(bin_size),chromosomes = chromosomes)

#load discriminated regions
dis_regions=read.table(dis_bedfile)

blacklist_haplotype=list()

x=bins_allele_m[[1]]
blacklist_haplotype[[allele_m_name]]=calc_black_list(x,dis_regions,allele_m_name)
x=bins_allele_p[[1]]
blacklist_haplotype[[allele_p_name]]=calc_black_list(x,dis_regions,allele_p_name)

##
save(blacklist_haplotype,file=paste0(outname,"_",allele_m_name,"_",allele_p_name,"_blacklist.Rdata"))
##

##########Export the blacklist filtered reads in control samples#######
haplo_control_S_G1_reads_newBL=list()
haplo_control_S_G1_reads_newBL[[allele_m_name]]=bed2GRanges(allele_m_bedfile,assembly=genome,remove.duplicate.reads = TRUE,min.mapq = 10,blacklist = blacklist_haplotype[[allele_m_name]])
haplo_control_S_G1_reads_newBL[[allele_p_name]]=bed2GRanges(allele_p_bedfile,assembly=genome,remove.duplicate.reads = TRUE,min.mapq = 10,blacklist = blacklist_haplotype[[allele_p_name]])

##
save(haplo_control_S_G1_reads_newBL,file=paste0(outname,"_",allele_m_name,"_",allele_p_name,"_control_G1_rmdup_reads_filterout_blacklist.Rdata"))
##

