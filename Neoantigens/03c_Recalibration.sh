#!/bin/bash

# Usage: 03_Recalibration.sh <output-folder> <sample_name>

# Input: BAM file, with duplicates tagged
# Output: BAM file, Quality scores for bases recalibrated

DIR_OUTPUT=$1
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$2


docker run --rm -v ${DIR_OUTPUT}:/Output \
        -v ${DIR_HG38}:/Hg38_dir \
        broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	ApplyBQSR \
	-R /Hg38_dir/Homo_sapiens_assembly38.fasta \
	-I /Output/${SAMPLE}.dedupped.bam \
	--bqsr-recal-file /Output/${SAMPLE}.perform_bqsr.table \
	-O /Output/${SAMPLE}.recal.bam

