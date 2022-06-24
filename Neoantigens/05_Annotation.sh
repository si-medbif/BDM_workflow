#!/bin/bash

# Usage: 05_Annotation.sh <output-folder> <sample_name>

# Input: VCF file, somatic mutations
# Output: VCF file, somatic mutations, annotated by VEP software package

# Required: vep docker image (see instructions in configure_vep.txt)

DIR_REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle
DIR_VCF=$1
SAMPLE=$2

docker run --rm -v ${DIR_REFERENCE}:/Reference \
	-v ${DIR_VCF}:/vcf \
	vep:v1 ./vep \
	--input_file /vcf/${SAMPLE}.m2.filtered.vcf.gz \
	--output_file /vcf/${SAMPLE}.vep.vcf \
	--format vcf --vcf --symbol --terms SO --tsl --hgvs \
	--fork 16 \
	--fasta /Reference/Homo_sapiens_assembly38.fasta \
	--offline --cache \
	--plugin Downstream --plugin Wildtype --force_overwrite

#./filter_annotated_vcf.py ${DIR_OUTPUT}/${SAMPLE}_m2_vep.vcf ${DIR_OUTPUT}/${SAMPLE}_m2_vep_filtered.vcf
