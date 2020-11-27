#!/bin/bash

BAMFOLDER=$1
dir_Hg38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$2


docker run --rm -v ${BAMFOLDER}:/Output \
	-v ${dir_Hg38}:/Hg38_dir \
	broadinstitute/gatk:4.1.9.0 gatk --java-options "-Xmx8G" \
	BaseRecalibrator \
	-R /Hg38_dir/Homo_sapiens_assembly38.fasta \
	--known-sites /Hg38_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites /Hg38_dir/dbsnp_146.hg38.vcf.gz \
	-I /Output/${SAMPLE}_dedupped.bam \
	-O /Output/${SAMPLE}_perform_bqsr.table

docker run --rm -v ${BAMFOLDER}:/Output \
        -v ${dir_Hg38}:/Hg38_dir \
        broadinstitute/gatk:4.1.9.0 gatk --java-options "-Xmx8G" \
	ApplyBQSR \
	-R /Hg38_dir/Homo_sapiens_assembly38.fasta \
	-I /Output/${SAMPLE}_dedupped.bam \
	--bqsr-recal-file /Output/${SAMPLE}_perform_bqsr.table \
	-O /Output/${SAMPLE}_recal.bam

