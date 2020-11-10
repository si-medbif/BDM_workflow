#!/bin/bash

# Running STAR alignment with parameters fit for using the results as input for STAR-fusion
# Parameter settings are based on recommendations for running STAR-fusion after alignment
# see: https://github.com/STAR-Fusion/STAR-Fusion/wiki

STARFOLDER=/home/harald/star
BAMFOLDER=$2
FASTQFOLDER=$1
HG38FOLDER=/gnome/WTS/GRCh38_v98
SAMPLE=$3

${STARFOLDER}/STAR \
	--runThreadN 4 \
	--genomeDir ${HG38FOLDER}/ \
        --readFilesIn ${FASTQFOLDER}/${SAMPLE}_1.fastq.gz ${FASTQFOLDER}/${SAMPLE}_2.fastq.gz \
	--outReadsUnmapped None \
        --readFilesCommand zcat \
        --outFileNamePrefix ${BAMFOLDER}/${SAMPLE}/$1_ \
        --outSAMtype BAM SortedByCoordinate \
	--outSAMstrandField intronMotif \
	--outSAMunmapped Within \
        --twopassMode Basic \
        --outSAMattributes NH HI AS nM XS \
	--chimSegmentMin 12 \
	--chimJunctionOverhangMin 8 \
	--chimOutJunctionFormat 1 \
	--alignSJDBoverhangMin 10 \
	--alignMatesGapMax 100000 \
	--alignIntronMax 100000 \
	--alignSJstitchMismatchNmax 5 -1 5 5 \
	--outSAMattrRGline ID:GRPundef \
	--chimMultimapScoreRange 3 \
	--chimScoreJunctionNonGTAG -4 \
	--chimMultimapNmax 20 \
	--chimNonchimScoreDropMin 10 \
	--peOverlapNbasesMin 12 \
	--peOverlapMMp 0.1 \
	--alignInsertionFlush Right \
	--alignSplicedMateMapLminOverLmate 0 \
	--alignSplicedMateMapLmin 30
