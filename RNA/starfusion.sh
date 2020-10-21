#! /bin/bash

set -euo pipefail

RESOURCE=/gnome/WTS/CTAT_Resource_Lib/GRCh38_v27_CTAT_lib_Feb092018/
FASTQ=$1
OUT=$2
FQSAMPLE=$3
SAMPLE=$4

docker run -v ${RESOURCE}:/resource \
	-v ${FASTQ}:/fastq \
	-v ${OUT}:/out \
	--rm trinityctat/starfusion \
	/usr/local/src/STAR-Fusion/STAR-Fusion \
	--left_fq /fastq/${FQSAMPLE}_R1.fastq.gz \
	--right_fq /fastq/${FQSAMPLE}_R2.fastq.gz \
	--genome_lib_dir /resource/ctat_genome_lib_build_dir \
	-O /out/starfusion/${SAMPLE} \
	--FusionInspector validate \
	--examine_coding_effect \
	--denovo_reconstruct
