#!/bin/bash

# Usage: 04_SomaticVariantCall.sh <output-folder> <GATK-resource-folder> <Tumor_sample> <Normal_sample>

# Input: BAM file, base quality scores are recalibrated.
# Output: VCF file, somatic variants, filtered.


DIR_VCF=$1
DIR_TUMORBAM=$2
DIR_NORMALBAM=$3
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
DIR_PON=/gnome/genome_database/pon_reference
TUMOR=$4
NORMAL=$5

docker run --rm -v ${DIR_VCF}:/vcf \
	-v ${DIR_NORMALBAM}:/normalbam \
	-v ${DIR_TUMORBAM}:/tumorbam \
        -v ${DIR_HG38}:/hg38 \
	-v ${DIR_PON}:/pon \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" Mutect2 \
        -R /hg38/Homo_sapiens_assembly38.fasta \
        -I /tumorbam/${TUMOR}.recal.bam \
        -tumor ${TUMOR} \
        -I /normalbam/${NORMAL}.recal.bam \
        -normal ${NORMAL} \
	-pon /pon/BTN.pon.vcf.gz \
	--f1r2-tar-gz /vcf/${TUMOR}.f1r2.tar.gz \
        --germline-resource /hg38/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /vcf/${TUMOR}.m2.vcf.gz

docker run --rm -v ${DIR_VCF}:/vcf \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" LearnReadOrientationModel \
	-I /vcf/${TUMOR}.f1r2.tar.gz \
	-O /vcf/${TUMOR}.read-orientation-model.tar.gz

docker run --rm -v ${DIR_HG38}:/hg38 \
	-v ${DIR_TUMORBAM}:/tumorbam \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /tumorbam/${TUMOR}.recal.bam \
	-V /hg38/small_exac_common_3.hg38.vcf.gz \
	-L /hg38/small_exac_common_3.hg38.vcf.gz \
	-O /tumorbam/${TUMOR}.getpileupsummaries.table

docker run --rm -v ${DIR_HG38}:/hg38 \
	-v ${DIR_NORMALBAM}:/normalbam \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" GetPileupSummaries \
	-I /normalbam/${NORMAL}.recal.bam \
	-V /hg38/small_exac_common_3.hg38.vcf.gz \
	-L /hg38/small_exac_common_3.hg38.vcf.gz \
	-O /normalbam/${NORMAL}.getpileupsummaries.table

docker run --rm -v ${DIR_TUMORBAM}:/tumorbam \
	-v ${DIR_NORMALBAM}:/normalbam \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" CalculateContamination \
	-I /tumorbam/${TUMOR}.getpileupsummaries.table \
	-tumor-segmentation /tumorbam/${TUMOR}.segments.table \
	-matched /normalbam/${NORMAL}.getpileupsummaries.table \
	-O /tumorbam/${TUMOR}.calculatecontamination.table

docker run --rm -v ${DIR_VCF}:/vcf \
	-v ${DIR_TUMORBAM}:/tumorbam \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-V /vcf/${TUMOR}.m2.vcf.gz \
	--contamination-table /tumorbam/${TUMOR}.calculatecontamination.table \
	--tumor-segmentation /tumorbam/${TUMOR}.segments.table \
	--ob-priors /vcf/${TUMOR}.read-orientation-model.tar.gz \
	-O /vcf/${TUMOR}.m2.filtered.vcf.gz
