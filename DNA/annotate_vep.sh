#!/bin/bash

# Cosmic annotation will be an ID on the form: COSV58737189
# For now, use this ID to look it up in the full database.

set -ueo pipefail

GATKhg38=/gnome/genome_database/gatk_bundle/hg38bundle
COSMIC=/gnome/genome_database/cosmic
DBNSFP=/gnome/genome_database/dbNSFP/dbNSFP4.0a
ENSEMBL=/gnome/genome_database/ensembl/vep_data
INPUT=$1
OUTPUT=$2
SAMPLE=$3

docker run --rm --name vep_${SAMPLE} \
	-v ${ENSEMBL}:/opt/vep/.vep \
	-v ${INPUT}:/input \
	-v ${OUTPUT}:/output \
	-v ${DBNSFP}:/dbnsfp \
	-v ${COSMIC}:/cosmic \
	-v ${GATKhg38}:/hg38 \
	ensemblorg/ensembl-vep \
	./vep --cache --offline --format vcf --vcf --force_overwrite \
	--dir_cache /opt/vep/.vep/ \
	--dir_plugins /opt/vep/.vep/Plugins/ \
	--input_file /output/${VCFFILE}.vcf.gz \
	--output_file /output/${VCFFILE}_vep.vcf \
	--stats_file /output/${VCFFILE}.summary.html \
	--plugin dbNSFP,/dbnsfp/dbNSFP4.0a.gz,ALL \
	--custom /cosmic/CosmicCodingMuts_chr_M.GRCh38.v94.vcf.gz,Cosmic,vcf,exact,0 \
	--assembly GRCh38 \
	--fork 8 \
	--species homo_sapiens \
	--cache_version 104 \
	--everything \
	--total_length

# Potential options, not used:
#	--filter_common : Remove variants with AF > 0.01
#	--per_gene : Only output the most severe consequence, arbitrary transcript.
#	--pick : Pick one line or block of consequence data per variant.
#	--flag_pick : As pick, but only adds a flag to the picked consequence.
