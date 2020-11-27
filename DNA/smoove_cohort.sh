#! /bin/bash

# Add additional bam files at the end following the established pattern

set -euo pipefail

PROJECT=$1
BAM=$2
OUT=$3
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE1=$4
SAMPLE2=$5

docker run --rm \
	-v ${BAM}:/bam \
	-v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	brentp/smoove smoove \
	call \
	--removepr \
	--name ${PROJECT} \
	--exclude /ref/exclude.cnvnator_100bp.GRCh38.20170403.bed \
	--fasta /ref/Homo_sapiens_assembly38.fasta \
	--outdir /out \
	--processes 16 \
	--genotype \
	/bam/${SAMPLE1}_recal.bam \
	/bam/${SAMPLE2}_recal.bam
