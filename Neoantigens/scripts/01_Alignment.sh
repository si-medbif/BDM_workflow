#!/bin/bash

# Usage: 01_alignment.sh <fastq-folder> <output-folder> <sample_name>

# Input: FASTQ-files on the form <sample_name>_R1.fastq.gz, <sample_name>_R2.fastq.gz
# Output: BAM file, sorted by coordinate
# Requires: samtools 
#           docker image biocontainers/bwa

dir_fastq=$1
dir_Output=$2
dir_Hg38=/gnome/genome_database/gatk_bundle/hg38bundle

mkdir -p ${dir_Output}/$3
mkdir -p ${dir_Output}/$3/{LOG,BAM,VCF}

docker run --rm -v ${dir_fastq}:/fastq \
	-v ${dir_Hg38}:/reference \
	biocontainers/bwa:v0.7.15_cv3 \
	bwa mem \
	-t 4 \
	-R "@RG\tID:$3\tSM:$3\tPL:Illumina\tLB:WES" \
	/reference/Homo_sapiens_assembly38.fasta.gz \
	/fastq/$3_R1.fastq.gz \
	/fastq/$3_R2.fastq.gz \
	| samtools sort -o ${dir_Output}/$3/BAM/$3_sorted.bam
