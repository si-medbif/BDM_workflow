#!/bin/bash

OUT=/gnome/chanitra_neoantigens/RNA2
BAM=/gnome/chanitra_neoantigens
REFERENCE=/gnome/genome_database/ensembl
DATABASE=/gnome/genome_database/gatk_bundle/hg38bundle
dir_GATK=/gnome/tutorial_datasets/gatk_somatic_variant_calling
dir_reference=/penguin/harald/neoantigens/chanitra_data

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
	-L /Reference/agilent_sureselect_V5_UTR_hg38_clean.bed \
        -pon /GATKtutorial/1000g_pon.hg38.vcf.gz \
        --germline-resource /GATKtutorial/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /out/$1_m2_raw.vcf

docker run --rm -v ${OUT}:/out \
	broadinstitute/gatk gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-V /out/$1_m2_raw.vcf \
	-O /out/$1_m2_filtered.vcf
