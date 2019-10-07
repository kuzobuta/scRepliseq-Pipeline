#!/usr/bin/env Rscript

#w200k s40k analysis AneuFinder#

library("pracma")
library("AneuFinder")

args = commandArgs(TRUE)
fragment_dir=args[1]
outdir=args[2]
blacklist=args[3]
genome_file=args[4]

dir.create(outdir,showWarnings=F)

window=200000
sliding=40000

files=list.files(fragment_dir,pattern=".Rdata",full.names=T)

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

for (fragment in files){
options(scipen=100)
name=basename(fragment)
name=sub("_mapq10_blacklist_fragment.Rdata", "", name)

load(fragment)
test=binReads(raw_reads,
              assembly=genome,
              blacklist=blacklist,
              bins=x_gr)

Total_reads=length(raw_reads)

test_read_counts=NULL
for (i in 1:(length(test)-1)){
  if(i==1){
    test_read_counts=data.frame(chr=seqnames(test[[i]]),starts=start(test[[i]])-1,ends=end(test[[i]]),c=test[[i]]$counts)
  }else{
  tmp=data.frame(chr=seqnames(test[[i]]),starts=start(test[[i]])-1,ends=end(test[[i]]),c=test[[i]]$counts)
  test_read_counts=rbind(test_read_counts,tmp)
  }
}

out=test_read_counts
a=out[order(out$chr,out$starts),]
a1=data.frame(chr=a[,1],s=round((a[,2]+a[,3])/2-sliding/2),e=round((a[,2]+a[,3])/2-sliding/2)+sliding,c=a$c/Total_reads)
a2=a1[!is.infinite(a1$c),]

write.table(a2,paste0(outdir,"/",name,"_w200ks40k_TagDensity.bedGraph"),col.names=F,row.names=F,quote=F,sep="\t")

options(scipen=0)
pdf(paste0(outdir,"/",name,"_w200ks40k_TagDensity.pdf"))
plot(density(a2$c[a2$c!=0]),main=paste0(name,"\nTotal reads:",Total_reads),xlab="Tag Density",xlim=c(0,0.0002))
dev.off()

}


