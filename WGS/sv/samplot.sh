#! /bin/bash

# Install samplot through conda: 
# $ conda install -c bioconda samplot 

TUMORBAM=$1
NORMALBAM=$2

TUMOR=$3
NORMAL=$4
CHROM=$5	# chr4
START=$6	# 160356811
END=$7		# 160358911
TYPE=$8		# DEL

samplot plot \
	-n ${TUMOR} ${NORMAL} \
	-b ${TUMORBAM}/${TUMOR}_recal.bam \
	  ${NORMALBAM}/${NORMAL}_recal.bam \
	-o ${TUMOR}.${CHROM}-${START}-${END}.png \
	-c ${CHROM} \
	-s ${START} \
	-e ${END} \
	-t ${TYPE}
