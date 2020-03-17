#!/bin/bash

STARFOLDER=/home/harald/star
BAMFOLDER=/gnome/tmp
FASTQFOLDER=/gnome/tmp
HG38FOLDER=/gnome/WTS/GRCh38_v98

${STARFOLDER}/STAR \
	--runThreadN 4 \
	--genomeDir ${HG38FOLDER}/ \
        --readFilesIn ${FASTQFOLDER}/$1_1.fastq.gz ${FASTQFOLDER}/$1_2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${BAMFOLDER}/$1_ \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --outSAMattributes NH HI AS nM XS 
