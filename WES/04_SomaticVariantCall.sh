#!/bin/bash

# Usage: 04_SomaticVariantCall.sh 
#           <output-folder> 
#           <GATK-resource-folder> 
#           <folder with bed file for WES regions> 
#           <Tumor_sample> 
#           <Normal_sample>

# Input: BAM file, base quality scores are recalibrated.
# Output: VCF file, somatic variants, filtered.


DIR_OUTPUT=$1
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
DIR_REFERENCE=$2
DIR_GATK=/gnome/tutorial_datasets/gatk_somatic_variant_calling/
TUMOR=$4
NORMAL=$5

docker run --rm -v ${DIR_OUTPUT}:/Output \
        -v ${DIR_HG38}:/Hg38_dir \
	-v ${DIR_GATK}:/GATKtutorial \
	-v ${DIR_REFERENCE}:/Reference \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" Mutect2 \
        -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
        -I /Output/${TUMOR}/BAM/${TUMOR}_recal.bam \
        -tumor ${TUMOR} \
        -I /Output/${NORMAL}/BAM/${NORMAL}_recal.bam \
        -normal ${NORMAL} \
	-L /Reference/agilent_sureselect_V5_UTR_hg38_clean.bed \
	-pon /GATKtutorial/1000g_pon.hg38.vcf.gz \
        --germline-resource /GATKtutorial/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /Output/${TUMOR}/VCF/${TUMOR}_m2.vcf.gz

docker run --rm -v ${DIR_GATK}:/data \
	-v ${DIR_OUTPUT}:/out \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /out/${TUMOR}/BAM/${TUMOR}_recal.bam \
	-V /data/small_exac_common_3.hg38.vcf.gz \
	-L /data/small_exac_common_3.hg38.vcf.gz \
	-O /out/${TUMOR}/BAM/tumor_getpileupsummaries.table

docker run --rm -v ${DIR_GATK}:/data \
	-v ${DIR_OUTPUT}:/out \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /out/${NORMAL}/BAM/${NORMAL}_recal.bam \
	-V /data/small_exac_common_3.hg38.vcf.gz \
	-L /data/small_exac_common_3.hg38.vcf.gz \
	-O /out/${TUMOR}/BAM/normal_getpileupsummaries.table

docker run --rm -v ${DIR_GATK}:/data \
	-v ${DIR_OUTPUT}:/out \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" CalculateContamination \
	-I /out/${TUMOR}/BAM/tumor_getpileupsummaries.table \
	-matched /out/${TUMOR}/BAM/normal_getpileupsummaries.table \
	-O /out/${TUMOR}/BAM/tumor_calculatecontamination.table

docker run --rm -v ${DIR_OUTPUT}:/out \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-V /out/${TUMOR}/VCF/${TUMOR}_m2.vcf.gz \
	--contamination-table /out/${TUMOR}/BAM/tumor_calculatecontamination.table \
	-O /out/${TUMOR}/VCF/${TUMOR}_m2_filtered.vcf.gz
