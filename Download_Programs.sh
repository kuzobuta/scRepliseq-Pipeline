#!/bin/sh

## SOFTWARE: bwa
## VERSION: 0.7.10
## TYPE: aligner
## SOURCE_URL: https://sourceforge.net/projects/bio-bwa/files/
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2
tar -xjf bwa-0.7.10.tar.bz2
cd bwa-0.7.10
make
cd ..
ln -s bwa-0.7.10 bwa

## SOFTWARE: samtools
## VERSION: 1.3.1
## TYPE: file format converter
## SOURCE_URL: https://github.com/samtools/samtools
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar -xjf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
make
cd ..
ln -s samtools-1.3.1 samtools
PATH=/usr/local/bin/samtools/:$PATH

## SOFTWARE: samstat
## VERSION: 1.5.1
## SOURCE_URL: https://sourceforge.net/projects/samstat
wget https://sourceforge.net/projects/samstat/files/samstat-1.5.1.tar.gz
tar -xzf samstat-1.5.1.tar.gz
cd samstat-1.5.1
PATH=/usr/local/bin/samtools/:$PATH
./configure
make
cd ..
ln -s samstat-1.5.1 samstat

## SOFTWARE: FastQC
wget https://github.com/s-andrews/FastQC/archive/v0.11.8.zip
unzip fastqc_v0.11.8.zip
chmod 755 FastQC/fastqc
ln -s FastQC/fastqc /usr/local/bin/fastqc
