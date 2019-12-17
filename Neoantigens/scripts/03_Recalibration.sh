#!/bin/bash

# Usage: 03_Recalibration.sh <output-folder> <sample_name>

# Input: BAM file, with duplicates tagged
# Output: BAM file, Quality scores for bases recalibrated

dir_Output=$1
dir_Hg38=/gnome/genome_database/gatk_bundle/hg38bundle
sample=$2


docker run --rm -v ${dir_Output}:/Output \
	-v ${dir_Hg38}:/Hg38_dir \
	broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	BaseRecalibrator \
	-R /Hg38_dir/Homo_sapiens_assembly38.fasta \
	--known-sites /Hg38_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites /Hg38_dir/dbsnp_146.hg38.vcf.gz \
	-I /Output/${sample}/BAM/${sample}_dedupped.bam \
	-O /Output/${sample}/BAM/${sample}_perform_bqsr.table

docker run --rm -v ${dir_Output}:/Output \
        -v ${dir_Hg38}:/Hg38_dir \
        broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	ApplyBQSR \
	-R /Hg38_dir/Homo_sapiens_assembly38.fasta \
	-I /Output/${sample}/BAM/${sample}_dedupped.bam \
	--bqsr-recal-file /Output/${sample}/BAM/${sample}_perform_bqsr.table \
	-O /Output/${sample}/BAM/${sample}_recal.bam

