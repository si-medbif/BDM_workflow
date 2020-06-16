#!/bin/bash

# Usage: 05_Annotation.sh <output-folder> <sample_name>

# Input: VCF file, somatic mutations
# Output: VCF file, somatic mutations, annotated by VEP software package

# Required: vep docker image (see instructions in configure_vep.txt)

DIR_REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle
VCFFOLDER=$1
SAMPLE=$2

docker run --rm -v ${DIR_REFERENCE}:/Reference \
	-v ${VCFFOLDER}:/Output \
	vep:v2 ./vep \
	--input_file /Output/${SAMPLE}_m2_filtered.vcf.gz \
	--output_file /Output/${SAMPLE}_m2_vep.vcf \
	--format vcf --vcf --symbol --terms SO --tsl --hgvs \
	--fasta /Reference/Homo_sapiens_assembly38.fasta \
	--offline --cache \
	--force_overwrite

