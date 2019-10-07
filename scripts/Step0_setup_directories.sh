#!/bin/bash

#scRepli-seq pipeline
#Step0 Setup directories

if [ $# -ne 1 ] ; then
    echo "Usage: bash Step0_setup_directories.sh [directory_name]"
    echo ""
    echo "Generate the directories for analyzing scRepliseq"
    echo "Example run:"
    echo "$ Step0_setup_directories.sh RPE1"
    echo "#Results of the directories
RPE1
├── Aneu_analysis
│   ├── G1_control
│   ├── HMM
│   │   ├── 1-somy
│   │   │   ├── Rdata
│   │   │   ├── binary
│   │   │   └── plot
│   │   ├── 2-somy
│   │   │   ├── Rdata
│   │   │   ├── binary
│   │   │   └── plot
│   │   └── repliscores
│   ├── Log2_Med
│   ├── MAD_score
│   ├── bins
│   └── fragment
├── bam
├── fastq
├── fastqc
├── logs
└── trim_fastq"
    echo ""
    exit 0
fi

result_dir=$1

mkdir -p ${result_dir}
mkdir -p ${result_dir}/fastq
mkdir -p ${result_dir}/fastqc
mkdir -p ${result_dir}/trim_fastq
mkdir -p ${result_dir}/bam
mkdir -p ${result_dir}/Aneu_analysis
mkdir -p ${result_dir}/Aneu_analysis/MAD_score
mkdir -p ${result_dir}/Aneu_analysis/Log2_Med
mkdir -p ${result_dir}/Aneu_analysis/bins
mkdir -p ${result_dir}/Aneu_analysis/fragment
mkdir -p ${result_dir}/Aneu_analysis/HMM
mkdir -p ${result_dir}/Aneu_analysis/HMM/1-somy
mkdir -p ${result_dir}/Aneu_analysis/HMM/2-somy
mkdir -p ${result_dir}/Aneu_analysis/HMM/1-somy/Rdata
mkdir -p ${result_dir}/Aneu_analysis/HMM/1-somy/binary
mkdir -p ${result_dir}/Aneu_analysis/HMM/1-somy/plot
mkdir -p ${result_dir}/Aneu_analysis/HMM/2-somy/Rdata
mkdir -p ${result_dir}/Aneu_analysis/HMM/2-somy/binary
mkdir -p ${result_dir}/Aneu_analysis/HMM/2-somy/plot
mkdir -p ${result_dir}/Aneu_analysis/HMM/repliscores
mkdir -p ${result_dir}/Aneu_analysis/G1_control
mkdir -p ${result_dir}/logs

