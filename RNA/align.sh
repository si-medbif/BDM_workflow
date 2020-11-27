#!/bin/bash

STARFOLDER=/home/harald/star
BAMFOLDER=$2
FASTQFOLDER=$1
HG38FOLDER=/gnome/WTS/GRCh38_v98
SAMPLE=$3

${STARFOLDER}/STAR \
	--runThreadN 4 \
	--genomeDir ${HG38FOLDER}/ \
        --readFilesIn ${FASTQFOLDER}/${SAMPLE}_1.fastq.gz ${FASTQFOLDER}/${SAMPLE}_2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${BAMFOLDER}/${SAMPLE}/$1_ \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outSAMattributes NH HI AS nM XS 
