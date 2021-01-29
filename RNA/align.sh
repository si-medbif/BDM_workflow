#!/bin/bash

STARFOLDER=/home/harald/star
BAMFOLDER=$2
FASTQFOLDER=$1
HG38FOLDER=/gnome/WTS/GRCh38_v98
SAMPLE=$3
LEFTFQ=$4
RIGHTFQ=$5

${STARFOLDER}/STAR \
	--runThreadN 4 \
	--genomeDir ${HG38FOLDER}/ \
        --readFilesIn ${FASTQFOLDER}/${LEFTFQ} ${FASTQFOLDER}/${RIGHTFQ} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${BAMFOLDER}/${SAMPLE}_ \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outSAMattributes NH HI AS nM XS 
