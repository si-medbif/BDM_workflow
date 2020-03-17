#!/bin/bash

OUT=/gnome/tmp
REFERENCE=/gnome/genome_database/ensembl
DATABASE=/gnome/genome_database/gatk_bundle/hg38bundle

#Recalibrate quality scores
docker run --rm -v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	-v ${DATABASE}:/db \
	broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	BaseRecalibrator \
	-R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	--known-sites /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--known-sites /db/dbsnp_146.hg38.vcf.gz \
	-I /out/$1_splitncigar.bam \
	-O /out/$1_perform_bqsr.table

docker run --rm -v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	ApplyBQSR \
	-R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	-I /out/$1_splitncigar.bam \
	--bqsr-recal-file /out/$1_perform_bqsr.table \
	-O /out/$1_recal.bam
