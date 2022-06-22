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
FQ_1=$4
FQ_2=$5

docker run --rm -v ${DIR_FASTQ}:/fastq \
	-v ${DIR_HG38}:/reference \
	biocontainers/bwa:v0.7.15_cv3 \
	bwa mem \
	-t 32 \
	-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina\tLB:WES" \
	/reference/Homo_sapiens_assembly38.fasta.gz \
	/fastq/${FQ_1} \
	/fastq/${FQ_2} \
	| samtools sort -o ${DIR_OUTPUT}/${SAMPLE}.sorted.bam
