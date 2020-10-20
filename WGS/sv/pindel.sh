#! /bin/bash

# Required: "pindelconfig.txt" must exist in the output folder
# Three columns: sample path, insert size, sample label
# --------------------------------------------------------------
# The path should be on the form /bam/{sample_name}.bam,
# "bam" refers to the internal folder in the docker container and 
# is mapped to th first input parameter to this script.

BAM_DIR=$1
Hg38_DIR=/gnome/genome_database/gatk_bundle/hg38bundle
OUT_DIR=$2
CONFIG=$3
PROJECT=$4
CHROM=$5

docker run --rm -v ${BAM_DIR}:/bam \
	-v ${Hg38_DIR}:/reference \
	-v ${OUT_DIR}:/out \
	shuangbroad/pindel:v0.2.5b8 pindel \
	-f /reference/Homo_sapiens_assembly38.fasta \
	-i /out/${CONFIG}.pindelconfig.txt \
	-c ${CHROM} \
	-o /out/${PROJECT} \
	-T 32

