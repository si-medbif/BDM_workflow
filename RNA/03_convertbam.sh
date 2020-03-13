#!/bin/bash


OUT=/gnome/chanitra_neoantigens/RNA2

# Add read group and sort SAM file
docker run --rm -v ${OUT}:/out \
	broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	AddOrReplaceReadGroups \
	-I /out/$1_Aligned.sortedByCoord.out.bam \
	-O /out/$1_rg_added_sorted.bam \
	-SO coordinate \
	-RGID $1 \
	-RGLB agilent \
	-RGPL illumina \
	-RGPU miseq \
	-RGSM $1

# Deduplicate
docker run --rm -v ${OUT}:/out \
	broadinstitute/picard:latest \
	MarkDuplicates \
	I=/out/$1_rg_added_sorted.bam \
	O=/out/$1_dedup.bam \
	CREATE_INDEX=true \
	M=/out/$1_dedup_output.metrics

