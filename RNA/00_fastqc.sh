#! /bin/bash

DIR_FASTQ=$1
SAMPLE=$2

docker run --rm -v ${DIR_FASTQ}:/fastq \
    gcatio/fastqc opt/FastQC/fastqc \
    -o /fastq \
    -t 2 \
    /fastq/${SAMPLE}_1.fq.gz \
    /fastq/${SAMPLE}_2.fq.gz
