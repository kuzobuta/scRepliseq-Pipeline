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
in_dirs_1=args[1]
in_dirs_2=args[2]
out_dir_m=args[3]
out_dir=args[4]
genome=args[5]
genome_file=args[6]
out_name=args[7]
max=args[8]
min=args[9]
sort_list_file=args[10]


if (is.na(sort_list_file)){
 sort_names=NULL
}else{
 tmp=read.table(sort_list_file,sep="\t",header=T)
 sort_names=tmp[,1]
}

#################For visualizing one group############

files=list.files(in_dirs_1)
files1=NULL
for (file in files){
 files1=c(files1,paste0(in_dirs_2,"/",file))
}

if (is.na(sort_list_file)){
 files1_srt=files1
}else{
 files1_srt=NULL
 for (n in sort_names){
   files1_srt=c(files1_srt,files1[grep(n,files1)])
 }
}

out_name=paste0(out_name,"_igv_session.xml")

IGV_XML_export(files1_srt,out_name,out_dir_m,out_dir,genome,max,min)
