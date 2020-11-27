#!/bin/bash

set -ueo pipefail

DBNSFP=/gnome/genome_database/dbNSFP/dbNSFP4.0a
ENSEMBL=/gnome/genome_database/ensembl/vep_data
OUTPUT=$1

docker run --rm -v ${ENSEMBL}:/opt/vep/.vep \
	-v ${OUTPUT}:/output \
	-v ${DBNSFP}:/dbnsfp \
	ensemblorg/ensembl-vep \
	./vep --cache --offline --format vcf --vcf --force_overwrite \
	--symbol --terms SO --tsl \
	--protein --sift b --polyphen b \
	--af --af_gnomad \
	--dir_cache /opt/vep/.vep/ \
	--dir_plugins /opt/vep/.vep/Plugins/ \
	--input_file /output/${VCFFILE}.vcf.gz \
	--output_file /output/${VCFFILE} \
	--plugin dbNSFP,/dbnsfp/dbNSFP4.0a.gz,ALL
