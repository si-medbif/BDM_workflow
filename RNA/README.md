##STAR

Alignment of pair-end reads are done by STAR (v2.6.1d).

Reference used is Ensembl (v98):
- Homo_sapiens.GRCh38.dna.primary_assembly.fa
- Homo_sapiens.GRCh38.98.gtf

Alignment is performed by doing a two pass run, using the junctions detected from the first pass to improve alignment in the second pass.

##FeatureCounts

Quantification of fragments are done by featureCounts (subread v1.6.3)

Used options:
- Counts fragments (-p)
- Both ends have to align (-B)
- Feature identifier is transcript ID (-g): transcript_id
- Use exons as the feature to count (-t exon)
- Minimum mapping quality 30 (-Q 30)
- Assign reads to all overlapping features (-O)

