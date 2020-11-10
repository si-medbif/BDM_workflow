#!/bin/bash

BAMFOLDER_TUMOR=$1
BAMFOLDER_NORMAL=$2
VCFFOLDER=$3
REFERENCE=/gnome/genome_database/ensembl
DIR_GATK=/gnome/tutorial_datasets/gatk_somatic_variant_calling
TUMOR=$4
NORMAL=$5

# --independent-mates was required for earlier gatk versions, needs testing
# Could potentially benefit from a region file

docker run --rm -v ${BAMFOLDER_TUMOR}:/bam_tumor \
	-v ${BAMFOLDER_NORMAL}:/bam_normal \
	-v ${VCFFOLDER}:/vcf_folder \
	-v ${DIR_GATK}:/GATKtutorial \
	-v ${REFERENCE}:/Reference \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" Mutect2 \
	-R /Reference/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	-I /bam_tumor/${TUMOR}_recal.bam \
	-tumor ${TUMOR} \
	-I /bam_normal/${NORMAL}_recal.bam \
	-normal ${NORMAL} \
	-pon /GATKtutorial/1000g_pon.hg38.vcf.gz \
	--germline-resource /GATKtutorial/af-only-gnomad.hg38.vcf.gz \
	--af-of-alleles-not-in-resource 0.0000025 \
	-O /vcf_folder/${TUMOR}_m2_raw.vcf.gz

docker run --rm -v ${dir_GATK}:/GATKtutorial \
	-v ${BAMFOLDER_TUMOR}:/bam_tumor \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /bam_tumor/${TUMOR}_recal.bam \
	-V /GATKtutorial/small_exac_common_3.hg38.vcf.gz \
	-L /GATKtutorial/small_exac_common_3.hg38.vcf.gz \
	-O /bam_tumor/${TUMOR}_getpileupsummaries.table

docker run --rm -v ${dir_GATK}:/GATKtutorial \
	-v ${BAMFOLDER_NORMAL}:/bam_normal \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /bam_normal/${NORMAL}_recal.bam \
	-V /GATKtutorial/small_exac_common_3.hg38.vcf.gz \
	-L /GATKtutorial/small_exac_common_3.hg38.vcf.gz \
	-O /bam_normal/${NORMAL}_getpileupsummaries.table

docker run --rm \
	-v ${BAMFOLDER_TUMOR}:/bam_tumor \
	-v ${BAMFOLDER_NORMAL}:/bam_normal \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" CalculateContamination \
	-I /bam_tumor/${TUMOR}_getpileupsummaries.table \
	-matched /bam_normal/${NORMAL}_getpileupsummaries.table \
	-O /bam_tumor/$1_tumor_calculatecontamination.table

docker run --rm -v ${VCFFOLDER}:/vcf_folder \
	-v ${BAMFOLDER_TUMOR}:/bam_tumor \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	-V /vcf_folder/${TUMOR}_m2_raw.vcf.gz \
	--contamination-table /bam_tumor/${TUMOR}_calculatecontamination.table \
	-O /vcf_folder/${TUMOR}_m2_filtered.vcf.gz
