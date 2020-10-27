#! /bin/bash

set -euo pipefail

CTAT_GENOME_LIB=/gnome/WTS/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir
FASTQ=$1
OUT=$2
SAMPLE=$3

docker run --rm -v ${FASTQ}:/fastq \
	-v ${OUT}:/out \
	-v ${CTAT_GENOME_LIB}:/ctat_genome_lib_dir \
	trinityctat/ctat_mutations \
	/usr/local/src/ctat-mutations/ctat_mutations \
	--left /fastq/${SAMPLE}_R1.fastq.gz \
	--right /fastq/${SAMPLE}_R2.fastq.gz \
	--out_dir /out/ctat_mutations_outdir \
	--threads 10 \
	--genome_lib_dir /ctat_genome_lib_dir
