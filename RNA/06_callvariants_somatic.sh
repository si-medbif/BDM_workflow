#!/bin/bash

OUT=/gnome/tmp
BAM=/gnome/tmp
REFERENCE=/gnome/genome_database/ensembl
DATABASE=/gnome/genome_database/gatk_bundle/hg38bundle
dir_GATK=/gnome/tutorial_datasets/gatk_somatic_variant_calling

docker run --rm -v ${OUT}:/out \
        -v ${DATABASE}:/db \
        -v ${dir_GATK}:/GATKtutorial \
        -v ${REFERENCE}:/ref \
	-v ${BAM}:/bam \
	-v ${dir_reference}:/Reference \
        broadinstitute/gatk gatk --java-options "-Xmx16g" Mutect2 \
        -R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
        -I /out/$1_recal.bam \
        -tumor $1 \
        -I /bam/$2/BAM/$2_recal.bam \
        -normal $2 \
        -pon /GATKtutorial/1000g_pon.hg38.vcf.gz \
        --germline-resource /GATKtutorial/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /out/$1_m2_raw.vcf

docker run --rm -v ${OUT}:/out \
	broadinstitute/gatk gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-V /out/$1_m2_raw.vcf \
	-O /out/$1_m2_filtered.vcf
