#!/usr/bin/env Rscript
path="/usr/local/bin"
#scRepli-seq pipeline
#Computation of log2 med score following mappability correction

library("pracma")
library("AneuFinder")
source(normalizePath(paste0(path,"/util/Aneufinder_Optional_script.R")))

args = commandArgs(TRUE)
fragment=args[1]
outdir=args[2]
ref_Rdata=args[3]
blacklist=args[4]
genome_file=args[5]

median_scale = function (datalist) {
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


window=200000
sliding=40000

name=basename(fragment)
name=sub("_mapq10_blacklist_fragment.Rdata", "", name)


##loading black list and genome Info##
genome_tmp <- read.table(genome_file,sep="\t") #
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

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

load(fragment)
test=binReads(raw_reads,
              assembly=genome,
              blacklist=blacklist,
              bins=x_gr)

ref_data = load(ref_Rdata)
control_100_S_G1_reads = get(ref_data)

test_map=list()
for (i in 1:(length(test)-1)){
test_map[[i]]=correctMappability(test[[i]],same.binsize = T,control_100_S_G1_reads,assembly=genome)
}

test_map_med_log2=NULL
for (i in 1:length(test_map)){
  test_map_med_log2=rbind(test_map_med_log2,median_scale(test_map[[i]]))
}

out=test_map_med_log2
a=out[order(out$chr,out$starts),]
a1=data.frame(chr=a[,1],s=round((a[,2]+a[,3])/2-sliding/2),e=round((a[,2]+a[,3])/2-sliding/2)+sliding,c=a$map_median_scale_log2)
a2=a1[!is.infinite(a1$c),]

write.table(a2,paste(outdir,"/",name,"_w200ks40k_map_count_median_log2.bedGraph",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
