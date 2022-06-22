#!/bin/bash

# Usage: 03_Recalibration.sh <output-folder> <sample_name>

# Input: BAM file, with duplicates tagged
# Output: BAM file, Quality scores for bases recalibrated

# https://gatk.broadinstitute.org/hc/en-us/community/posts/360056396392-BaseRecalibrator-takes-forever-to-run-Any-suggestions-
# https://gatk.broadinstitute.org/hc/en-us/articles/360035889551-When-should-I-restrict-my-analysis-to-specific-intervals-
# https://gatk.broadinstitute.org/hc/en-us/articles/360056968612-GatherBQSRReports

DIR_OUTPUT=$1
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$2
CHROM=$3


docker run --rm -v ${DIR_OUTPUT}:/Output \
	-v ${DIR_HG38}:/Hg38_dir \
	broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	BaseRecalibrator \
	-R /Hg38_dir/Homo_sapiens_assembly38.fasta \
	--known-sites /Hg38_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites /Hg38_dir/dbsnp_146.hg38.vcf.gz \
	-I /Output/${SAMPLE}.dedupped.bam \
	-O /Output/${SAMPLE}.perform_bqsr.${CHROM}.table \
	-L ${CHROM}


