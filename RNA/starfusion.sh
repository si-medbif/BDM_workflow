#! /bin/bash

set -euo pipefail

RESOURCE=/gnome/WTS/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play
FASTQ=$1
OUT=$2
FQSAMPLE1=$3
FQSAMPLE2=$4
SAMPLE=$5

docker run -v ${RESOURCE}:/resource \
	-v ${FASTQ}:/fastq \
	-v ${OUT}:/out \
	--rm trinityctat/starfusion \
	/usr/local/src/STAR-Fusion/STAR-Fusion \
	--left_fq /fastq/${FQSAMPLE1} \
	--right_fq /fastq/${FQSAMPLE2} \
	--genome_lib_dir /resource/ctat_genome_lib_build_dir \
	-O /out/${SAMPLE} \
	--FusionInspector validate \
	--examine_coding_effect \
	--denovo_reconstruct
