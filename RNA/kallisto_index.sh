#!/bin/bash

REFERENCE=/gnome/genome_database/ensembl

docker run --rm -v ${REFERENCE}:/ref \
	zlskidmore/kallisto:latest kallisto index \
	-i /ref/Homo_sapiens.GRCh38.cdna.all.fa.kallisto.indexed \
	/ref/Homo_sapiens.GRCh38.cdna.all.fa.gz 
