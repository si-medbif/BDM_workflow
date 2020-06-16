#!/bin/bash

# Usage: 01_alignment.sh <fastq-folder> <output-folder> <sample_name>

# Input: FASTQ-files on the form <sample_name>_R1.fastq.gz, <sample_name>_R2.fastq.gz
# Output: BAM file, sorted by coordinate
# Requires: samtools 
#           docker image biocontainers/bwa

DIR_FASTQ=$1
DIR_OUTPUT=$2
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$3

mkdir -p ${dir_Output}/${SAMPLE}
mkdir -p ${dir_Output}/${SAMPLE}/{BAM,VCF}

docker run --rm -v ${DIR_FASTQ}:/fastq \
	-v ${DIR_HG38}:/reference \
	biocontainers/bwa:v0.7.15_cv3 \
	bwa mem \
	-t 4 \
	-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina\tLB:WES" \
	/reference/Homo_sapiens_assembly38.fasta.gz \
	/fastq/${SAMPLE}_R1.fastq.gz \
	/fastq/${SAMPLE}_R2.fastq.gz \
	| samtools sort -o ${DIR_OUTPUT}/${SAMPLE}/BAM/${SAMPLE}_sorted.bam
