#!/bin/bash

# Usage: 04_SomaticVariantCall.sh <output-folder> <GATK-resource-folder> <Tumor_sample> <Normal_sample>

# WGS: 0.19m 0.01m 13.24m 6.39m

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
	I=/vcf/${TUMOR}.m2.chr1.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr2.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr3.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr4.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr5.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr6.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr7.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr8.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr9.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr10.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr11.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr12.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr13.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr14.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr15.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr16.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr17.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr18.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr19.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr20.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr21.vcf.gz \
	I=/vcf/${TUMOR}.m2.chr22.vcf.gz \
	O=/vcf/${TUMOR}.m2.vcf.gz


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

#all_f1r2_input=`for chromosome in {1..22}; do
#        printf -- "-I ${TUMOR}.f1r2.chr${chromosome}tar.gz "; done`

docker run --rm -v ${DIR_VCF}:/vcf \
	broadinstitute/gatk:4.1.5.0 gatk \
	--java-options "-Xmx8g" LearnReadOrientationModel \
	-I /vcf/${TUMOR}.f1r2.chr1.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr2.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr3.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr4.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr5.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr6.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr7.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr8.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr9.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr10.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr11.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr12.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr13.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr14.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr15.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr16.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr17.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr18.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr19.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr20.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr21.tar.gz \
	-I /vcf/${TUMOR}.f1r2.chr22.tar.gz \
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
	--ob-priors /vcf/${TUMOR}.read-orientation-model.tar.gz \
	-O /vcf/${TUMOR}.m2.filtered.vcf.gz

mkdir ${DIR_VCF}/temp
mv ${DIR_VCF}/${TUMOR}.m2.chr*.vcf.gz.stats ${DIR_VCF}/temp
mv ${DIR_VCF}/${TUMOR}.m2.chr*.vcf.gz* ${DIR_VCF}/temp
mv ${DIR_VCF}/${TUMOR}.f1r2.chr*.tar.gz ${DIR_VCF}/temp

