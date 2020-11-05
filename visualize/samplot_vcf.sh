#! /bin/bash

# Install samplot through conda: 
# $ conda install -c bioconda samplot 

TUMORBAM=$1
NORMALBAM=$2
OUT=$3
VCF=$4

samplot vcf \
	--filter "SVTYPE == 'DEL' & SU >= 1" \
	--filter "SVTYPE == 'DUP' & SU >= 1" \
	--filter "SVTYPE == 'INV' & SU >= 1" \
	--vcf ${VCF} \
	-d ${OUT}/ \
	-O png \
	-T /gnome/genome_database/ensembl/Homo_sapiens.GRCh38.101.sort.gff3.gz \
	-A /gnome/genome_database/ucsc/rmsk.bed.gz \
	-b ${TUMORBAM} \
	${NORMALBAM} \
	> samplot_commands.sh
