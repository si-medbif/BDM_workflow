#! /bin/bash

DIR_FASTQ=$1
SAMPLE=$2

docker run --rm -v ${DIR_FASTQ}:/fastq \
    gcatio/fastqc opt/FastQC/fastqc \
    -o /fastq \
    -t 2 \
    /fastq/${SAMPLE}_R1.fastq.gz \
    /fastq/${SAMPLE}_R2.fastq.gz
