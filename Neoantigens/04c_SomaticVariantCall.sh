#!/bin/bash

# Usage: 04_SomaticVariantCall.sh <output-folder> <GATK-resource-folder> <Tumor_sample> <Normal_sample>

# Input: BAM file, base quality scores are recalibrated.
# Output: VCF file, somatic variants, filtered.

set -e

DIR_VCF=$1
DIR_TUMORBAM=$2
DIR_NORMALBAM=$3
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
DIR_PON=/gnome/genome_database/pon_reference
TUMOR=$4
NORMAL=$5

docker run --rm -v ${DIR_VCF}:/vcf \
	broadinstitute/picard:latest \
	MergeVcfs \
	I=/vcf/TBi_D01.m2.chr1.vcf.gz \
	I=/vcf/TBi_D01.m2.chr2.vcf.gz \
	I=/vcf/TBi_D01.m2.chr3.vcf.gz \
	I=/vcf/TBi_D01.m2.chr4.vcf.gz \
	I=/vcf/TBi_D01.m2.chr5.vcf.gz \
	I=/vcf/TBi_D01.m2.chr6.vcf.gz \
	I=/vcf/TBi_D01.m2.chr7.vcf.gz \
	I=/vcf/TBi_D01.m2.chr8.vcf.gz \
	I=/vcf/TBi_D01.m2.chr9.vcf.gz \
	I=/vcf/TBi_D01.m2.chr10.vcf.gz \
	I=/vcf/TBi_D01.m2.chr11.vcf.gz \
	I=/vcf/TBi_D01.m2.chr12.vcf.gz \
	I=/vcf/TBi_D01.m2.chr13.vcf.gz \
	I=/vcf/TBi_D01.m2.chr14.vcf.gz \
	I=/vcf/TBi_D01.m2.chr15.vcf.gz \
	I=/vcf/TBi_D01.m2.chr16.vcf.gz \
	I=/vcf/TBi_D01.m2.chr17.vcf.gz \
	I=/vcf/TBi_D01.m2.chr18.vcf.gz \
	I=/vcf/TBi_D01.m2.chr19.vcf.gz \
	I=/vcf/TBi_D01.m2.chr20.vcf.gz \
	I=/vcf/TBi_D01.m2.chr21.vcf.gz \
	I=/vcf/TBi_D01.m2.chr22.vcf.gz \
	O=/vcf/TBi_D01.m2.vcf.gz


docker run --rm -v ${DIR_VCF}:/vcf \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" MergeMutectStats \
	-stats /vcf/${TUMOR}.m2.chr1.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr2.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr3.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr4.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr5.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr6.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr7.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr8.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr9.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr10.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr11.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr12.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr13.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr14.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr15.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr16.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr17.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr18.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr19.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr20.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr21.vcf.gz.stats \
	-stats /vcf/${TUMOR}.m2.chr22.vcf.gz.stats \
	-O /vcf/${TUMOR}.m2.vcf.gz.stats

docker run --rm -v ${DIR_VCF}:/vcf \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" LearnReadOrientationModel \
	-I /vcf/${TUMOR}.f1r2.tar.gz \
	-O /vcf/${TUMOR}.read-orientation-model.tar.gz

docker run --rm -v ${DIR_VCF}:/vcf \
	-v ${DIR_TUMORBAM}:/tumorbam \
	-v ${DIR_HG38}:/hg38 \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-R /hg38/Homo_sapiens_assembly38.fasta \
	-V /vcf/${TUMOR}.m2.vcf.gz \
	--contamination-table /tumorbam/${TUMOR}.calculatecontamination.table \
	--filtering-stats /vcf/${TUMOR}.merged.stats \
	-O /vcf/${TUMOR}.m2.filtered.vcf.gz

mkdir temp
mv ${DIR_VCF}/${TUMOR}.m2.chr*.vcf.gz.stats /vcf/temp
mv ${DIR_VCF}/${TUMOR}.m2.chr*.vcf.gz* /vcf/temp
mv ${DIR_VCF}/${TUMOR}.f1r2.chr*.tar.gz /vcf/temp

