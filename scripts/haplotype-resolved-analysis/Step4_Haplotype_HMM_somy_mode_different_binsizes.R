#!/usr/bin/env Rscript

#v190403 Haplotype HMM analysis Step4
#
#
#$ Rscript Step4_Haplotype_HMM_somy_mode_different_binsizes.R <genome_file> \
#  <R_haplotype_file> \
#  <allele_m_name> <allele_p_name> \
#  <bin_size> \
#  <outdir> \
#  <Ref_file> \
#  <black_file> \
#  <somy>

library(AneuFinder)
options(scipen=100)
source("Aneufinder_Optional_script.R")
args = commandArgs(TRUE)

###Genome Info###
genome_tmp <- read.table(args[1],sep="\t") #UCSC_mm9.woYwR.fa.fai
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

R_file=args[2] #From step1 as haplo_bins_reads : ${outname}_haplotype_bins.Rdata
allele_m_name=args[3] #cba
allele_p_name=args[4] #msm
binsize=args[5] #400000
out_dir=args[6] # 
Ref_file=args[7] #From step 3 as haplo_control_S_G1_reads_newBL : ${outname}_${allele_m_name}_${allele_p_name}_control_G1_rmdup_reads_filterout_blacklist.Rdata
black_file=args[8] #blacklist_haplotype : ${outname}_${allele_m_name}_${allele_p_name}_blacklist.Rdata
somy=args[9] #1-somy or 2-somy

name=basename(R_file)
name=sub("_haplotype_bins.Rdata","",name)

out_dir_Rdata=paste(out_dir,"/","Rdata",sep="")
out_dir_binary=paste(out_dir,"/","binary",sep="")
out_dir_plot=paste(out_dir,"/","plot",sep="")

dir.create(out_dir,showWarnings = F)
dir.create(out_dir_Rdata,showWarnings = F)
dir.create(out_dir_plot,showWarnings = F)
dir.create(out_dir_binary,showWarnings = F)
dir.create(paste(out_dir_binary,"/",allele_m_name[1],sep=""),showWarnings = F)
dir.create(paste(out_dir_binary,"/",allele_p_name[1],sep=""),showWarnings = F)


load(Ref_file) 
load(black_file) #blacklist_haplotype

haplo_MAP_2HMM=list()
load(R_file)

allele_m=haplo_bins_reads[[allele_m_name]][[binsize]]
allele_p=haplo_bins_reads[[allele_p_name]][[binsize]]

ov_mat=findOverlaps(allele_m,blacklist_haplotype[[allele_m_name]])
ov_pat=findOverlaps(allele_p,blacklist_haplotype[[allele_p_name]])

mat_filter=allele_m
pat_filter=allele_p

#Blacklist regions to 0 counts
mat_filter$counts[queryHits(ov_mat)]=0
mat_filter$pcounts[queryHits(ov_mat)]=0
mat_filter$mcounts[queryHits(ov_mat)]=0

pat_filter$counts[queryHits(ov_pat)]=0
pat_filter$pcounts[queryHits(ov_pat)]=0
pat_filter$mcounts[queryHits(ov_pat)]=0

options(warn=1)
mat.MAP = correctMappability(mat_filter,
                             same.binsize=TRUE,
                             haplo_control_S_G1_reads_newBL[[allele_m_name]],
                             assembly=genome)
pat.MAP = correctMappability(pat_filter,
                             same.binsize=TRUE,
                             haplo_control_S_G1_reads_newBL[[allele_p_name]],
                             assembly=genome)

options(warn=2)
haplo_MAP_2HMM[[allele_m_name]]=findCNVs(mat.MAP,
                                         most.frequent.state = somy,
                                         method = 'HMM',num.threads=10,
                                         eps = 0.01,
                                         states=c("zero-inflation","0-somy","1-somy","2-somy"),
                                         max.iter = 5000,
                                         count.cutoff.quantile = 1)
haplo_MAP_2HMM[[allele_p_name]]=findCNVs(pat.MAP,
                                         most.frequent.state = somy,
                                         method = 'HMM',num.threads=10,
                                         eps = 0.01,
                                         states=c("zero-inflation","0-somy","1-somy","2-somy"),
                                         max.iter = 5000,
                                         count.cutoff.quantile = 1)


haplo_MAP_2HMM[[allele_m_name]]$ID=paste(name,allele_m_name,sep="_")
haplo_MAP_2HMM[[allele_p_name]]$ID=paste(name,allele_p_name,sep="_")

##save as Rdata
size=as.numeric(binsize)/1000
save(haplo_MAP_2HMM,file=paste(out_dir_Rdata,"/",name,"_",size,"k_2HMM_eps0.01_rmdup_newblack_",somy,".Rdata",sep=""))

##Export to bedGraph##
#allele_maternal#
out=data.frame(chr=seqnames(haplo_MAP_2HMM[[allele_m_name]]$bins),
               starts=start(haplo_MAP_2HMM[[allele_m_name]]$bins)-1,
               ends=end(haplo_MAP_2HMM[[allele_m_name]]$bins),
               copy.numbers=elementMetadata(haplo_MAP_2HMM[[allele_m_name]]$bins)$copy.number)
binary_out=out
binary_out$copy.numbers[out$copy.numbers==1]=-1
binary_out$copy.numbers[out$copy.numbers==2]=1
write.table(binary_out,paste(out_dir_binary,"/",allele_m_name,"/",name,"_",allele_m_name,"_",size,"k_100S_MAP_HMM2_eps0.01.binary_",somy,".bedGraph",sep=""),sep="\t",col.names=F,row.names=F,quote=F)

#allele_paternal#
out=data.frame(chr=seqnames(haplo_MAP_2HMM[[allele_p_name]]$bins),
               starts=start(haplo_MAP_2HMM[[allele_p_name]]$bins)-1,
               ends=end(haplo_MAP_2HMM[[allele_p_name]]$bins),
               copy.numbers=elementMetadata(haplo_MAP_2HMM[[allele_p_name]]$bins)$copy.number)
binary_out=out
binary_out$copy.numbers[out$copy.numbers==1]=-1
binary_out$copy.numbers[out$copy.numbers==2]=1
write.table(binary_out,paste(out_dir_binary,"/",allele_p_name,"/",name,"_",allele_p_name,"_",size,"k_100S_MAP_HMM2_eps0.01.binary_",somy,".bedGraph",sep=""),sep="\t",col.names=F,row.names=F,quote=F)



##export the plots##
#allele_maternal#
pdf(paste(out_dir_plot,"/",name,"_",allele_m_name,"_",size,"k_2HMM_eps0.01_",somy,".pdf",sep=""))
print(plot(haplo_MAP_2HMM$cba,type=2))
dev.off()

#allele_paternal#
pdf(paste(out_dir_plot,"/",name,"_",allele_p_name,"_",size,"k_2HMM_eps0.01_",somy,".pdf",sep=""))
print(plot(haplo_MAP_2HMM$msm,type=2))
dev.off()


