#!/bin/bash

# Usage: 04_SomaticVariantCall.sh 
#           <output-folder> 
#           <GATK-resource-folder> 
#           <folder with bed file for WES regions> 
#           <Tumor_sample> 
#           <Normal_sample>

# Input: BAM file, base quality scores are recalibrated.
# Output: VCF file, somatic variants, filtered.


BAMFOLDER_TUMOR=$1
BAMFOLDER_NORMAL=$2
VCFFOLDER=$3
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
BEDFOLDER=$4
DIR_GATK=/gnome/tutorial_datasets/gatk_somatic_variant_calling/
TUMOR=$5
NORMAL=$6

docker run --rm -v ${BAMFOLDER_TUMOR}:/bam_tumor \
	-v ${BAMFOLDER_NORMAL}:/bam_normal \
	-v ${VCFFOLDER}:/vcf_folder \
        -v ${DIR_HG38}:/Hg38_dir \
	-v ${DIR_GATK}:/GATKtutorial \
	-v ${BEDFOLDER}:/Reference \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" Mutect2 \
        -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
        -I /bam_tumor/${TUMOR}_recal.bam \
        -tumor ${TUMOR} \
        -I /bam_normal/${NORMAL}_recal.bam \
        -normal ${NORMAL} \
	-L /Reference/agilent_sureselect_V5_UTR_hg38_clean.bed \
	-pon /GATKtutorial/1000g_pon.hg38.vcf.gz \
        --germline-resource /GATKtutorial/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /vcf_folder/${TUMOR}_m2.vcf.gz

docker run --rm -v ${DIR_GATK}:/data \
	-v ${BAMFOLDER_TUMOR}:/out \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /out/${TUMOR}_recal.bam \
	-V /data/small_exac_common_3.hg38.vcf.gz \
	-L /data/small_exac_common_3.hg38.vcf.gz \
	-O /out/${TUMOR}_getpileupsummaries.table

docker run --rm -v ${DIR_GATK}:/data \
	-v ${BAMFOLDER_NORMAL}:/out \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /out/${NORMAL}_recal.bam \
	-V /data/small_exac_common_3.hg38.vcf.gz \
	-L /data/small_exac_common_3.hg38.vcf.gz \
	-O /out/${NORMAL}_getpileupsummaries.table

docker run --rm -v ${DIR_GATK}:/data \
	-v ${BAMFOLDER_TUMOR}:/bam_tumor \
	-v ${BAMFOLDER_NORMAL}:/bam_normal \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" CalculateContamination \
	-I /bam_tumor/${TUMOR}_getpileupsummaries.table \
	-matched /bam_normal/${NORMAL}_getpileupsummaries.table \
	-O /bam_tumor/${TUMOR}_calculatecontamination.table

docker run --rm -v ${VCFFOLDER}:/out \
	-v ${BAMFOLDER_TUMOR}:/bam_folder \
	-v ${DIR_HG38}:/Hg38_dir \
	broadinstitute/gatk:4.1.9.0 gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-R /Hg38_dir/Homo_sapiens_assembly38.fasta \
	-V /out/${TUMOR}_m2.vcf.gz \
	--contamination-table /bam_folder/${TUMOR}_calculatecontamination.table \
	-O /out/${TUMOR}_m2_filtered.vcf.gz
