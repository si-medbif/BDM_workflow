#!/bin/bash

set -e

STARFOLDER=/home/harald/star
ENSEMBL=/gnome/genome_database/ensembl
OUTREF=/gnome/WTS

## NB! STAR does not create any directories by itself!

mkdir -p ${OUTREF}/GRCh38_v98

# Create genome index files for STAR aligner
${STARFOLDER}/STAR \
	--runThreadN 20 \
	--runMode genomeGenerate \
	--genomeDir ${OUTREF}/GRCh38_v98/ \
	--genomeFastaFiles ${ENSEMBL}/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	--sjdbGTFfile ${ENSEMBL}/Homo_sapiens.GRCh38.98.convert.gtf \
	--sjdbOverhang 99

# Create reference index for GATK
samtools faidx ${ENSEMBL}/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa

# Create reference dictionary for GATK
docker run --rm -v ${ENSEMBL}:/ensembl \
	broadinstitute/picard:latest \
	CreateSequenceDictionary \
	R=/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	O=/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.convert.dict
