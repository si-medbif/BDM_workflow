#! /bin/bash

DIR_FASTQ=$1
FQ_LEFT=$2
FQ_RIGHT=$3


docker run --rm -v ${DIR_FASTQ}:/fastq \
    gcatio/fastqc opt/FastQC/fastqc \
    -o /fastq \
    -t 2 \
    /fastq/${FQ_LEFT} \
    /fastq/${FQ_RIGHT}
