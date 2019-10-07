#!/usr/bin/env Rscript
path="/usr/local/bin"

#############################################################
#Generation of IGV session for visualizing multiple binarized scRepli-seq data
#In IGV, it will show the binarized scRepli-seq data in heatmap mode (blue: replicated (+1), yellow: unreplicated (-1))
#Output file name should be name as "XXXX.xml"
#If you have a lot of samples (>50), we prefer to convert the bedGraph file to bigWig format by "bedGraphToBigWig" command 
#

options(scipen = 100)
source(normalizePath(paste0(path,"/util/IGV_Session_XML_one_two.R")))


args = commandArgs(TRUE)
in_dirs_1_m=args[1] #log2 med
in_dirs_1=args[2]
in_dirs_2_m=args[3] #binary
in_dirs_2=args[4]
out_dir_m=args[5]
out_dir=args[6]
genome=args[7]
genome_file=args[8]
out_name=args[9]
sort_list_file=args[10]

if (is.na(sort_list_file)){
 sort_names=NULL
}else{
 tmp=read.table(sort_list_file,sep="\t",header=T)
 sort_names=tmp[,1]
}


#################For visualizing two groups############
#To separate the data sets , you can generate the blank bed file based on your genome
genome_tmp <- read.table(genome_file,sep="\t") #UCSC_mm9.woYwR.fa.fai
blank_bed_file=cbind(as.character(genome_tmp$V1),rep(0,dim(genome_tmp)[1]),genome_tmp$V2)


write.table(blank_bed_file,
            paste0(out_dir_m,"/",genome,"_blank.bed"),col.names=F,row.names=F,quote=F,sep="\t")

blankfile=paste0(out_dir,"/",genome,"_blank.bed")

files=list.files(in_dirs_1_m)
files1=NULL
for (file in files){
 files1=c(files1,paste0(in_dirs_1,"/",file))
}

if (is.na(sort_list_file)){
 files1_srt=files1
}else{
 files1_srt=NULL
 for (n in sort_names){
   files1_srt=c(files1_srt,files1[grep(n,files1)])
 }
}

max1=0.75 #RT scores
min1=-0.75 #RT scores

files=list.files(in_dirs_2_m)
files2=NULL
for (file in files){
 files2=c(files2,paste0(in_dirs_2,"/",file))
}

if (is.na(sort_list_file)){
 files2_srt=files2
}else{
 files2_srt=NULL
 for (n in sort_names){
   files2_srt=c(files2_srt,files2[grep(n,files2)])
 }
}


max2=1 #binary data
min2=-1 #binary data

out_name=paste0(out_name,"_igv_session.xml")

IGV_XML_export_two_set(files1_srt,files2_srt,out_name,out_dir_m,out_dir,genome,blankfile,max1,min1,max2,min2)
