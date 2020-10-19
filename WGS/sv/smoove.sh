#! /bin/bash

set -euo pipefail

PROJECT=$1
BAM=$2
OUT=$3
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle


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
	/bam/BB-T0006-DNA_recal.bam \
	/bam/BB-T0012-DNA_recal.bam \
	/bam/BB-T0006-DNA_recal.bam
#	/bam/*.bam
