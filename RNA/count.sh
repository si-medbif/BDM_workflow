#!/bin/bash
# Counts fragments (-p) from pairs where both ends align (-B)
# Specify feature identifier gene_id/gene_name (-g gene_id/gene_name/transcript_id)
# Use exons as the feature to count (-t exon)
# (Potential) Perform read counting at feature level (-f) (e.g. exons)
# Minimum mapping quality (-Q <int>)
# Assign reads to all overlapping features (-O)
# Other useful options:
# -J  Count number of reads supporting each exon-exon junction. Results are saved to a file named '<output_file>.jcounts' 

FEATURECOUNTS=/home/harald/subread-1.6.3-Linux-x86_64/bin
BAMFOLDER=$1
GTF=/gnome/genome_database/ensembl/Homo_sapiens.GRCh38.98.convert.gtf
OUTFOLDER=$2
SAMPLE=$3

${FEATURECOUNTS}/featureCounts \
        -p \
        -B \
        -t exon \
        -g gene_id \
        -T 1 \
	-Q 30 \
	-O \
        -a ${GTF} \
        -o ${OUTFOLDER}/${SAMPLE}.readcounts.featurecounts.txt \
        ${BAMFOLDER}/${SAMPLE}_Aligned.sortedByCoord.out.bam 
