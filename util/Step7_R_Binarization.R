#!/usr/bin/env Rscript
path="/usr/local/bin"
#Binarization by HMM

library(AneuFinder)
source(normalizePath(paste0(path,"/util/Aneufinder_Optional_script.R")))
options(scipen = 100)

args = commandArgs(TRUE)
binfile=args[1]
out_dir=args[2]
Ref_file=args[3]
genome_file=args[4]
binsize=args[5]
somy=args[6]

##loading black list and genome Info##
genome_tmp <- read.table(genome_file,sep="\t") #
genome=data.frame(UCSC_seqlevel=genome_tmp$V1,UCSC_seqlength=genome_tmp$V2)
chromosomes=as.character(genome$UCSC_seqlevel)

dir.create(paste0(out_dir,"/binary"),showWarnings = FALSE)
dir.create(paste0(out_dir,"/Rdata"),showWarnings = FALSE)
dir.create(paste0(out_dir,"/plot"),showWarnings = FALSE)

name=basename(binfile)
name=sub("_mapq10_blacklist_bin.Rdata", "", name)

tmp=load(Ref_file)
control_S_G1_reads=get(tmp)

load(binfile)

binsize=as.character(binsize)
options(warn=1)
binned.MAP = correctMappability(bins_reads[[binsize]],
                             same.binsize=TRUE,
                             control_S_G1_reads,
                             assembly=genome)

if (is.na(somy)){
 somy="2-somy"
}

#options(warn=2)
#Try 3 run and export the file name
HMM_err=0
#1st try
t = try(findCNVs(binned.MAP,
                 ID=name,
                 method = 'HMM',num.threads=10,
                 eps = 0.01,
                 most.frequent.state = somy,
                 states=c("zero-inflation","0-somy","1-somy","2-somy"),
                 max.iter = 3000)
        ,silent=T)
if(class(t) == "try-error"){
  #2nd try
  HMM_err=HMM_err + 1
  t = try(findCNVs(binned.MAP,
                   ID=name,
                   method = 'HMM',num.threads=10,
                   eps = 0.01,
                   most.frequent.state = somy,
                   states=c("zero-inflation","0-somy","1-somy","2-somy"),
                   max.iter = 3000)
          ,silent=T)
  if(class(t) == "try-error"){
    #3rd try
    t = try(findCNVs(binned.MAP,
                     ID=name,
                     method = 'HMM',num.threads=10,
                     eps = 0.01,
                     most.frequent.state = somy,
                     states=c("zero-inflation","0-somy","1-somy","2-somy"),
                     max.iter = 3000)
            ,silent=T)
    if(class(t) == "try-error"){
      cat(name,"\t",out_dir,"\tDid not work by 3 runs","\n",file="Error_file_for_2HMM.txt",append=TRUE)
      stop("Did not work by 3 runs!!")
    }
  }
}

binned.MAP.HMM.2=t
#binned.MAP.HMM.2$ID=c(name)

ksize=as.numeric(binsize)/1000

#Save Rdata (HMM results)
save(binned.MAP.HMM.2,file=paste0(out_dir,"/Rdata/",name,"_",ksize,"k_100S_MAP_2HMM_eps0.01_",somy,".Rdata"))

#Plot density
pdf(file=paste0(out_dir,"/plot/",name,"_",ksize,"k_100S_MAP_2HMM_eps0.01_",somy,".pdf"))
print(plot(binned.MAP.HMM.2,type=2))
dev.off()

#Export binary bedGraph
out=data.frame(chr=seqnames(binned.MAP.HMM.2[[2]]),starts=start(binned.MAP.HMM.2[[2]])-1,ends=end(binned.MAP.HMM.2[[2]]),copy.numbers=elementMetadata(binned.MAP.HMM.2[[2]])$copy.number)
binary_out=out
binary_out$copy.numbers[out$copy.numbers==1]=-1
binary_out$copy.numbers[out$copy.numbers==2]=1
write.table(binary_out,paste0(out_dir,"/binary/",name,"_",ksize,"k_100S_MAP_2HMM2_eps0.01.binary_",somy,".bedGraph"),sep="\t",col.names=F,row.names=F,quote=F)



