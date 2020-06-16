#!/bin/bash

dir_Output=$1
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle
CHROM=$2

docker run --rm -v ${dir_Output}:/Output \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.5.0 gatk --java-options "-Xmx64g" \
	GenotypeGVCFs \
	-R /ref/Homo_sapiens_assembly38.fasta \
	-V gendb:///Output/project_db.${CHROM} \
	-O /Output/project.${CHROM}.vcf.gz \
	--tmp-dir=/Output/tmp
