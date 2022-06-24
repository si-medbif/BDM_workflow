#!/bin/bash

# Usage: 04_SomaticVariantCall.sh <output-folder> <GATK-resource-folder> <Tumor_sample> <Normal_sample>

# Input: BAM file, base quality scores are recalibrated.
# Output: VCF file, somatic variants, filtered.


DIR_VCF=$1
TUMOR=$2

docker run --rm -v ${DIR_VCF}:/vcf \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" SelectVariants \
	-V /vcf/${TUMOR}.vep.vcf \
	-O /vcf/${TUMOR}.vep.filtered.vcf \
	--exclude-filtered
