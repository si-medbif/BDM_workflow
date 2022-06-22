#!/bin/bash

# Usage: 04_SomaticVariantCall.sh <output-folder> <GATK-resource-folder> <Tumor_sample> <Normal_sample>

# WES: 66.78m (chr1)
# WGS: 233.65m (chr1)

# Input: BAM file, base quality scores are recalibrated.
# Output: VCF file, somatic variants, filtered.


DIR_VCF=$1
DIR_TUMORBAM=$2
DIR_NORMALBAM=$3
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
DIR_PON=/gnome/genome_database/pon_reference
TUMOR=$4
NORMAL=$5
CHROM=$6

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
	-L ${CHROM} \
	-pon /pon/BTN.pon.vcf.gz \
	--f1r2-tar-gz /vcf/${TUMOR}.f1r2.${CHROM}.tar.gz \
        --germline-resource /hg38/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /vcf/${TUMOR}.m2.${CHROM}.vcf.gz

