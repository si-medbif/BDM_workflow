#!/bin/bash

BAM_TUMOR=$1
BAM_NORMAL=$2
REFERENCE=/gnome/genome_database/ensembl
DATABASE=/gnome/genome_database/gatk_bundle/hg38bundle
dir_GATK=/gnome/tutorial_datasets/gatk_somatic_variant_calling
TUMOR=$3
NORMAL=$4

docker run --rm -v ${BAM_TUMOR}:/out \
        -v ${DATABASE}:/db \
        -v ${dir_GATK}:/GATKtutorial \
        -v ${REFERENCE}:/ref \
	-v ${BAM_NORMAL}:/bam \
	-v ${dir_reference}:/Reference \
        broadinstitute/gatk gatk --java-options "-Xmx16g" Mutect2 \
        -R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
        -I /out/${TUMOR}/${TUMOR}_recal.bam \
        -tumor ${TUMOR} \
        -I /bam/${NORMAL}/BAM/${NORMAL}_recal.bam \
        -normal ${NORMAL} \
        -pon /GATKtutorial/1000g_pon.hg38.vcf.gz \
        --germline-resource /GATKtutorial/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /out/${TUMOR}_m2_raw.vcf

docker run --rm -v ${BAM_TUMOR}:/out \
	broadinstitute/gatk gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-V /out/${TUMOR}/${TUMOR}_m2_raw.vcf \
	-O /out/${TUMOR}/${TUMOR}_m2_filtered.vcf
