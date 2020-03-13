#!/bin/bash
# Counts fragments (-p) from pairs where both ends align (-B)
# Specify feature identifier gene_id/gene_name (-g gene_id/gene_name/transcript_id)
# Use exons as the feature to count (-t exon)
# (Potential) Perform read counting at feature level (-f) (e.g. exons)
# Minimum mapping quality (-Q <int>)
# Assign reads to all overlapping features (-O)
FEATURECOUNTS=/home/harald/subread-1.6.3-Linux-x86_64/bin
BAMFOLDER=/gnome/chanitra_neoantigens/RNA2
GTF=/gnome/genome_database/ensembl/Homo_sapiens.GRCh38.98.convert.gtf
OUTFOLDER=/gnome/chanitra_neoantigens/RNA2

${FEATURECOUNTS}/featureCounts \
        -p \
        -B \
        -t exon \
        -g transcript_id \
        -T 1 \
	-Q 30 \
	-O \
        -a ${GTF} \
        -o ${OUTFOLDER}/readcounts.featurecounts.txt \
        ${BAMFOLDER}/$1_Aligned.sortedByCoord.out.bam 
