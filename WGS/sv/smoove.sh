#! /bin/bash

set -euo pipefail

PROJECT=$1
TUMORBAM=$2
NORMALBAM=$3
TUMOR=$4
NORMAL=$5
OUT=$6
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle


docker run --rm \
	-v ${TUMORBAM}:/tumor \
	-v ${NORMALBAM}:/normal \
	-v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	brentp/smoove smoove \
	call \
	--removepr \
	--name ${PROJECT} \
	--exclude /ref/exclude.cnvnator_100bp.GRCh38.20170403.bed \
	--fasta /ref/Homo_sapiens_assembly38.fasta \
	--outdir /out \
	--processes 8 \
	--genotype \
	/tumor/${TUMOR}_recal.bam \
	/normal/${NORMAL}_recal.bam
