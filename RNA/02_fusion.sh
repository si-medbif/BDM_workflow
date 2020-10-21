#! /bin/bash

# From https://github.com/STAR-Fusion/STAR-Fusion/wiki

set -euo pipefail

RESOURCE=/gnome/WTS/CTAT_Resource_Lib/GRCh38_v27_CTAT_lib_Feb092018/
ALIGNFOLDER=$1
OUT=$2
SAMPLE=$3

docker run -v ${RESOURCE}:/resource \
	-v ${ALIGNFOLDER}:/align \
	-v ${OUT}:/out \
	--rm trinityctat/starfusion \
	/usr/local/src/STAR-Fusion/STAR-Fusion \
	--genome_lib_dir /resource/ctat_genome_lib_build_dir \
	-J /align/Chimeric.out.junction \
	-O /out/${SAMPLE}
