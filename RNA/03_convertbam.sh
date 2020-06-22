#!/bin/bash


OUT=$1
SAMPLE=$2

# Add read group and sort SAM file
docker run --rm -v ${OUT}:/out \
	broadinstitute/gatk:4.1.7.0 gatk --java-options "-Xmx8G" \
	AddOrReplaceReadGroups \
	-I /out/${SAMPLE}_Aligned.sortedByCoord.out.bam \
	-O /out/${SAMPLE}_rg_added_sorted.bam \
	-SO coordinate \
	-RGID ${SAMPLE} \
	-RGLB agilent \
	-RGPL illumina \
	-RGPU miseq \
	-RGSM ${SAMPLE}

# Deduplicate
docker run --rm -v ${OUT}:/out \
	broadinstitute/picard:latest \
	MarkDuplicates \
	I=/out/${SAMPLE}_rg_added_sorted.bam \
	O=/out/${SAMPLE}_dedup.bam \
	CREATE_INDEX=true \
	M=/out/${SAMPLE}_dedup_output.metrics

