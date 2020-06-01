#!/bin/bash

OUT=$1
REFERENCE=/gnome/genome_database/ensembl
DATABASE=/gnome/genome_database/gatk_bundle/hg38bundle
dir_GATK=/gnome/tutorial_datasets/gatk_somatic_variant_calling
SAMPLE=$2

docker run --rm -v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	-v ${DATABASE}:/db \
	broadinstitute/gatk gatk --java-options "-Xmx8G" \
	HaplotypeCaller \
	-R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	-I /out/${SAMPLE}_recal.bam \
	--dont-use-soft-clipped-bases \
	-stand-call-conf 20.0 \
	-O /out/${SAMPLE}_hc_raw.vcf 

#Filter variants
# Hard filtering due to lack of good (any) references for VQSR
docker run --rm -v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	-v ${DATABASE}:/db \
	broadinstitute/gatk gatk --java-options "-Xmx8G" \
	VariantFiltration \
	-R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	-V /out/${SAMPLE}_hc_raw.vcf \
	-O /out/${SAMPLE}_hc_filtered.vcf \
	-window 35 \
	-cluster 3 \
	--filter-name FS \
	-filter "FS > 30.0" \
	--filter-name QD \
	-filter "QD < 2.0"

