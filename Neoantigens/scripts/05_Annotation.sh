#!/bin/bash

# Usage: 05_Annotation.sh <output-folder> <sample_name>

# Input: VCF file, somatic mutations
# Output: VCF file, somatic mutations, annotated by VEP software package

# Required: vep docker image (see instructions in configure_vep.txt)

dir_Reference=/gnome/genome_database/gatk_bundle/hg38bundle
dir_Output=$1
sample=$2

docker run --rm -v ${dir_Reference}:/Reference \
	-v ${dir_Output}:/Output \
	vep:v2 ./vep \
	--input_file /Output/${sample}/VCF/${sample}_m2_filtered.vcf.gz \
	--output_file /Output/${sample}/VCF/${sample}_m2_vep.vcf \
	--format vcf --vcf --symbol --terms SO --tsl --hgvs \
	--fasta /Reference/Homo_sapiens_assembly38.fasta \
	--offline --cache \
	--plugin Downstream --plugin Wildtype --force_overwrite
