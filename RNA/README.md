## STAR

Alignment of pair-end reads are done by STAR (v2.6.1d).

Reference used is Ensembl (v98):
- Homo_sapiens.GRCh38.dna.primary_assembly.fa
- Homo_sapiens.GRCh38.98.gtf

Alignment is performed by doing a two pass run, using the junctions detected from the first pass to improve alignment in the second pass.

## FeatureCounts

Quantification of fragments are done by featureCounts (subread v1.6.3)

Used options:
- Counts fragments (-p)
- Both ends have to align (-B)
- Feature identifier is transcript ID (-g): transcript_id
- Use exons as the feature to count (-t exon)
- Minimum mapping quality 30 (-Q 30)
- Assign reads to all overlapping features (-O)

## Suggested workflow for measuring expression levels:

1. Create a project folder.
2. Within the project folder, create a FASTQ folder, a BAM folder and a COUNT folder.
3. Copy all the FASTQ-files into the fastq-folder.
4. Do qualitycontrol on fastq files:
	1. ./00_fastqc.sh FASTQFOLDER SAMPLE
	2. ./create_fastqc_report.py FASTQFOLDER SAMPLE1 SAMPLE2 ...
5. Create BAM files:
	1. ./01_align.sh FASTQFOLDER BAMFOLDER SAMPLE
	2. ./create_bam_report.py BAMFOLDER SAMPLE1 SAMPLE2 ...
6. Count expression levels:
	1. ./02_count.sh BAMFOLDER OUTFOLDER SAMPLE
	2. ./count2tpm.py OUTFOLDER
7. Combine the created figures into a final report.
	1. ./create_final_report.py

## Suggested workflow for detecting variants

1. Follow expression levels workflow until creating of BAM files
2. Prepare BAM files for variant detection
	1. ./03_convertbam.sh BAMFOLDER SAMPLE
	2. ./04_splitncigar.sh BAMFOLDER SAMPLE
	3. ./05_recalibrate.sh BAMFOLDER SAMPLE
3. Call variants
	1. ./06_callvariants.sh BAMFOLDER SAMPLE
	-- or --
	1. ./06_callvariants_somatic.sh BAMFOLDER_TUMOR BAMFOLDER_NORMAL SAMPLE_TUMOR SAMPLE_NORMAL

