#!/bin/bash

# Usage: 01_alignment.sh <sample_name>

# Input: FASTQ-files on the form <sample_name>_R1.fastq.gz, <sample_name>_R2.fastq.gz
# Requires: samtools v.1.7
#           docker image biocontainers/bwa

FASTQFOLDER=$1
BAMFOLDER=$2
dir_Hg38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$3
LEFTFQ=$4
RIGHTFQ=$5

docker run --rm -v ${FASTQFOLDER}:/fastq \
	-v ${dir_Hg38}:/reference \
	biocontainers/bwa:v0.7.15_cv3 \
	bwa mem \
	-t 16 \
	-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina\tLB:WGS" \
	/reference/Homo_sapiens_assembly38.fasta.gz \
	/fastq/${LEFTFQ} \
	/fastq/${RIGHTFQ} \
	| samtools sort -o ${BAMFOLDER}/${SAMPLE}_sorted.bam
