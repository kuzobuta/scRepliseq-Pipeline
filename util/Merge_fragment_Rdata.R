#!/usr/bin/env Rscript

#Merge the fragment Rdata

library("optparse")
library("AneuFinder")

option_list = list(
	make_option(c("-o", "--out"), type="character", default="Merged_fragment.Rdata", 
                 help="output file name [default= %default]", metavar="character")) 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser,positional_arguments=TRUE); 

files=opt$args

tmp_raw_reads=NULL
i=0
for (file in files){
 i=i+1
 load(file)
 if(i==1){
  tmp_raw_reads=raw_reads
 }else{
  tmp_raw_reads=c(tmp_raw_reads,raw_reads)
 }
}

out_file=opt$options$out
raw_reads=tmp_raw_reads
save(raw_reads,file=out_file)

