# Workflow

## FASTQC
```
00_Fastqc.sh $FASTQ ${SAMPLE1}_R1.fastq.gz ${SAMPLE1}_R2.fastq.gz
```
```
docker run --rm -v ${DIR_FASTQ}:/fastq gcatio/fastqc opt/FastQC/fastqc -o /fastq -t 2 /fastq/${FQ_LEFT} /fastq/${FQ_RIGHT}
```

## RNAseq analysis

run_rnaseq.sh
```
nextflow run nf-core/rnaseq \
        -profile docker \
        --input input.csv \
        --outdir ${PWD} \
        --igenomes_base ${DATABASES}/ngi-igenomes \
        --genome GRCh38 \
        -resume
```

## HLA-estimation:
run_optitype.sh
```
nextflow run nf-core/hlatyping -profile docker --input '${FASTQ}/*_R{1,2}.fastq.gz'
```
run_arcashla.sh
```
./docker_arcashla.sh ${SAMPLE}
```
```
docker run \
        --rm -v ${BAMFOLDER}:/bam \
        -v ${FASTQ}:/fastq \
        arcashla:v1 \
        arcasHLA extract \
        -t 16 \
        --paired \
        --o /fastq \
        /bam/${SAMPLE}.markdup.sorted.bam

docker run \
        --rm -v ${FASTQ}:/fastq \
        -v ${OUT}:/out \
        arcashla:v1 \
        arcasHLA genotype \
        -t 16 \
        --outdir /out \
        /fastq/${SAMPLE}.markdup.sorted.extracted.1.fq.gz \
        /fastq/${SAMPLE}.markdup.sorted.extracted.2.fq.gz
```

run_hlahd.sh
```
docker_hlahd.sh \
	${TUMOR} \
	${FASTQ_FOLDER} \
	${TUMOR}_R1.fastq.gz \
	${TUMOR}_R2.fastq.gz \
	${OUT}
```
```
docker run --rm -v ${DIR_FASTQ}:/fastq \
	-v ${OUT}:/out \
	jianhung/hla-hd hlahd.sh \
	-t 16 \
	-f /home/hlahd.1.2.1/freq_data \
	fastq/${FASTQ1} \
	fastq/${FASTQ2} \
	/home/hlahd.1.2.1/HLA_gene.split.txt \
	/home/hlahd.1.2.1/dictionary \
	 ${SAMPLE} \
	 /out
```

## Align reads
run_alignment.sh
```
01_Alignment.sh ${DIR_FASTQ} ${DIR_OUTPUT} ${SAMPLE1} ${FQ_1} ${FQ_2}
```
```
docker run --rm -v ${DIR_FASTQ}:/fastq \
        -v ${DIR_HG38}:/reference \
        biocontainers/bwa:v0.7.15_cv3 \
        bwa mem \
        -t 32 \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina\tLB:WES" \
        /reference/Homo_sapiens_assembly38.fasta.gz \
        /fastq/${FQ_1} \
        /fastq/${FQ_2} \
        | samtools sort -o ${DIR_OUTPUT}/${SAMPLE}.sorted.bam
```

## Mark duplicates
run_markduplicate.sh
```
02_MarkDuplicate.sh ${DIR_OUTPUT} ${SAMPLE1} &
```
```
docker run --rm -v ${DIR_OUTPUT}:/Output \
        broadinstitute/picard:latest \
        MarkDuplicates \
        CREATE_INDEX=true \
        I=/Output/${SAMPLE}.sorted.bam \
        O=/Output/${SAMPLE}.dedupped.bam \
        M=/Output/${SAMPLE}.dedup_output.metrics
```

## Recalibrate samples
run_recalibration1.sh
```
for i in {1..22}; do
        CHROM=chr$i
        03a_Recalibration.sh ${BAMFOLDER} ${SAMPLE} ${CHROM} &
        echo "/Output/${SAMPLE}.perform_bqsr.${CHROM}.table" >> ${BAMFOLDER}/${SAMPLE}.perform_bqsr.list
done
```
```
docker run --rm -v ${DIR_OUTPUT}:/Output \
        -v ${DIR_HG38}:/Hg38_dir \
        broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
        BaseRecalibrator \
        -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
        --known-sites /Hg38_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --known-sites /Hg38_dir/dbsnp_146.hg38.vcf.gz \
        -I /Output/${SAMPLE}.dedupped.bam \
        -O /Output/${SAMPLE}.perform_bqsr.${CHROM}.table \
        -L ${CHROM}
```

run_recalibration2.sh
```
03b_Recalibration.sh ${BAMFOLDER} ${SAMPLE1}
```
```
docker run --rm -v ${DIR_OUTPUT}:/Output \
        -v ${DIR_HG38}:/Hg38_dir \
        broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
        BaseRecalibrator \
        -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
        --known-sites /Hg38_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --known-sites /Hg38_dir/dbsnp_146.hg38.vcf.gz \
        -I /Output/${SAMPLE}.dedupped.bam \
        -O /Output/${SAMPLE}.perform_bqsr.${CHROM}.table \
        -L ${CHROM}
```
run_recalibration3.sh
```
03c_Recalibration.sh ${BAMFOLDER} ${SAMPLE1} &
```
```
docker run --rm -v ${DIR_OUTPUT}:/Output \
        -v ${DIR_HG38}:/Hg38_dir \
        broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
        GatherBQSRReports \
        -I /Output/${SAMPLE}.perform_bqsr.list \
        -O /Output/${SAMPLE}.perform_bqsr.table
		
docker run --rm -v ${DIR_OUTPUT}:/Output \
        -v ${DIR_HG38}:/Hg38_dir \
        broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
        ApplyBQSR \
        -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
        -I /Output/${SAMPLE}.dedupped.bam \
        --bqsr-recal-file /Output/${SAMPLE}.perform_bqsr.table \
        -O /Output/${SAMPLE}.recal.bam
```

## Call somatic variants
run_somaticvariants1.sh
```
for i in {1..22}; do
        CHROM=chr${i}
        04a_SomaticVariantCall.sh ${DIR_VCF} ${DIR_TUMORBAM} ${DIR_NORMALBAM} ${TUMOR} ${NORMAL} ${CHROM} &
done
```
```
docker run --rm -v ${DIR_VCF}:/vcf \
        -v ${DIR_NORMALBAM}:/normalbam \
        -v ${DIR_TUMORBAM}:/tumorbam \
        -v ${DIR_HG38}:/hg38 \
        -v ${DIR_PON}:/pon \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" Mutect2 \
        -R /hg38/Homo_sapiens_assembly38.fasta \
        -I /tumorbam/${TUMOR}.recal.bam \
        -tumor ${TUMOR} \
        -I /normalbam/${NORMAL}.recal.bam \
        -normal ${NORMAL} \
        -L ${CHROM} \
        -pon /pon/BTN.pon.vcf.gz \
        --f1r2-tar-gz /vcf/${TUMOR}.f1r2.${CHROM}.tar.gz \
        --germline-resource /hg38/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /vcf/${TUMOR}.m2.${CHROM}.vcf.gz
```
run_somaticvariants2.sh
```
04b_SomaticVariantCall.sh ${DIR_VCF} ${DIR_TUMORBAM} ${DIR_NORMALBAM} ${TUMOR} ${NORMAL} &
```
```
docker run --rm -v ${DIR_HG38}:/hg38 \
        -v ${DIR_TUMORBAM}:/tumorbam \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" GetPileupSummaries \
        -I /tumorbam/${TUMOR}.recal.bam \
        -V /hg38/small_exac_common_3.hg38.vcf.gz \
        -L /hg38/small_exac_common_3.hg38.vcf.gz \
        -O /tumorbam/${TUMOR}.getpileupsummaries.table

docker run --rm -v ${DIR_HG38}:/hg38 \
        -v ${DIR_NORMALBAM}:/normalbam \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" GetPileupSummaries \
        -I /normalbam/${NORMAL}.recal.bam \
        -V /hg38/small_exac_common_3.hg38.vcf.gz \
        -L /hg38/small_exac_common_3.hg38.vcf.gz \
        -O /normalbam/${NORMAL}.getpileupsummaries.table

docker run --rm -v ${DIR_TUMORBAM}:/tumorbam \
        -v ${DIR_NORMALBAM}:/normalbam \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" CalculateContamination \
        -I /tumorbam/${TUMOR}.getpileupsummaries.table \
        -tumor-segmentation /tumorbam/${TUMOR}.segments.table \
        -matched /normalbam/${NORMAL}.getpileupsummaries.table \
        -O /tumorbam/${TUMOR}.calculatecontamination.table
```

run_somaticvariants3.sh
```
04c_SomaticVariantCall.sh ${DIR_VCF} ${DIR_TUMORBAM} ${DIR_NORMALBAM} ${TUMOR} ${NORMAL}
```		
```
docker run --rm -v ${DIR_VCF}:/vcf \
        broadinstitute/picard:latest \
        MergeVcfs \
        I=/vcf/${TUMOR}.m2.chr1.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr2.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr3.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr4.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr5.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr6.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr7.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr8.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr9.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr10.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr11.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr12.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr13.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr14.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr15.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr16.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr17.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr18.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr19.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr20.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr21.vcf.gz \
        I=/vcf/${TUMOR}.m2.chr22.vcf.gz \
        O=/vcf/${TUMOR}.m2.vcf.gz
		
docker run --rm -v ${DIR_VCF}:/vcf \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" MergeMutectStats \
        -stats /vcf/${TUMOR}.m2.chr1.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr2.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr3.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr4.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr5.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr6.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr7.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr8.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr9.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr10.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr11.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr12.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr13.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr14.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr15.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr16.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr17.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr18.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr19.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr20.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr21.vcf.gz.stats \
        -stats /vcf/${TUMOR}.m2.chr22.vcf.gz.stats \
        -O /vcf/${TUMOR}.m2.vcf.gz.stats
		
docker run --rm -v ${DIR_VCF}:/vcf \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" LearnReadOrientationModel \
        -I /vcf/${TUMOR}.f1r2.chr1.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr2.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr3.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr4.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr5.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr6.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr7.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr8.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr9.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr10.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr11.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr12.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr13.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr14.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr15.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr16.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr17.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr18.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr19.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr20.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr21.tar.gz \
        -I /vcf/${TUMOR}.f1r2.chr22.tar.gz \
        -O /vcf/${TUMOR}.read-orientation-model.tar.gz
		
docker run --rm -v ${DIR_VCF}:/vcf \
        -v ${DIR_TUMORBAM}:/tumorbam \
        -v ${DIR_HG38}:/hg38 \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" FilterMutectCalls \
        -R /hg38/Homo_sapiens_assembly38.fasta \
        -V /vcf/${TUMOR}.m2.vcf.gz \
        --contamination-table /tumorbam/${TUMOR}.calculatecontamination.table \
        --filtering-stats /vcf/${TUMOR}.merged.stats \
        --ob-priors /vcf/${TUMOR}.read-orientation-model.tar.gz \
        -O /vcf/${TUMOR}.m2.filtered.vcf.gz

mkdir -p ${DIR_VCF}/temp
mv ${DIR_VCF}/${TUMOR}.m2.chr*.vcf.gz.stats ${DIR_VCF}/temp
mv ${DIR_VCF}/${TUMOR}.m2.chr*.vcf.gz* ${DIR_VCF}/temp
mv ${DIR_VCF}/${TUMOR}.f1r2.chr*.tar.gz ${DIR_VCF}/temp
```

## Annotate variants
run_annotate.sh
```
05_Annotation.sh ${DIR_VCF} ${TUMOR} &
```
```
docker run --rm -v ${DIR_REFERENCE}:/Reference \
        -v ${DIR_VCF}:/vcf \
        vep:v2 ./vep \
        --input_file /vcf/${SAMPLE}.m2.filtered.vcf.gz \
        --output_file /vcf/${SAMPLE}.vep.vcf \
        --format vcf --vcf --symbol --terms SO --tsl --hgvs \
        --fork 16 \
        --fasta /Reference/Homo_sapiens_assembly38.fasta \
        --offline --cache \
        --plugin Frameshift --plugin Wildtype --force_overwrite
```

## Add coverage
run_addcoverage.sh
```
07_AddCoverage.sh ${DIR_VCF} ${DIR_RNABAM} ${DIR_TPM} ${RNASAMPLE} ${TUMOR} &
```
```
docker run --rm -v ${DIR_VCF}:/vcf \
        -v ${DIR_HG38}:/Hg38_dir \
        -v ${DIR_RNABAM}:/rnabam \
        mgibio/bam_readcount_helper-cwl:1.0.0 \
        /usr/bin/python /usr/bin/bam_readcount_helper.py \
        /vcf/${TUMOR}.vep.vcf \
        ${TUMOR} \
        /Hg38_dir/Homo_sapiens_assembly38.fasta \
        /rnabam/${RNASAMPLE}.markdup.sorted.bam \
        /vcf/

# Add gene expression data to the vcf files
docker run --rm -v ${DIR_VCF}:/vcf \
        -v ${DIR_TPM}:/tpm \
        griffithlab/vatools \
        vcf-expression-annotator \
        /vcf/${TUMOR}.vep.vcf \
        /tpm/${RNASAMPLE}.tpm.tsv \
        custom \
        gene \
        --id-column gene_id \
        --expression-column ${RNASAMPLE} \
        -s ${TUMOR}

# Add RNA read counts to the vcf files
docker run --rm -v ${DIR_VCF}:/vcf \
        vatools:v1 \
        vcf-readcount-annotator \
        /vcf/${TUMOR}.vep.gx.vcf \
        /vcf/${TUMOR}_bam_readcount_snv.tsv \
        RNA \
        -s ${TUMOR} \
        -t snv \
        -o /vcf/${TUMOR}_1.vcf

docker run --rm -v ${DIR_VCF}:/vcf \
        vatools:v1 \
        vcf-readcount-annotator \
        /vcf/${TUMOR}_1.vcf \
        /vcf/${TUMOR}_bam_readcount_indel.tsv \
        RNA \
        -s ${TUMOR} \
        -t indel \
        -o /vcf/${TUMOR}.vep.gx.cov.vcf
```

## Select variants
run_selectvariants.sh
```
06_SelectVariants.sh ${DIR_VCF} ${IN_VCF} ${OUT_VCF}
```
```
docker run --rm -v ${DIR_VCF}:/vcf \
        broadinstitute/gatk:4.1.5.0 gatk \
        --java-options "-Xmx8g" SelectVariants \
        -V /vcf/${IN_VCF} \
        -O /vcf/${OUT_VCF} \
        --exclude-filtered
```

## Predict neoantigens
run_predictepitopes.sh
```
HLA=HLA-A*02:03,HLA-B*40:06,HLA-C*07:06,DPA1*01:03,DPB1*04:02,DQA1*01:02,DQB1*02:02,DRA*01:01,DRB1*07:01
mkdir -p ${dir_PVAC}/${TUMOR}
08_EpitopePrediction_new.sh ${dir_VCF} ${dir_PVAC} ${TUMOR} ${NORMAL} ${HLA}
```
```
docker run --rm -v ${dir_VCF}:/vcf \
       -v ${dir_PVAC}/${TUMOR}:/pvac \
       griffithlab/pvactools:3.0.1 \
        pvacseq run \
       -t 16 \
       -e1 8,9,10,11,12 \
       --iedb-install-directory /opt/iedb \
       -a sample_name \
       --normal-sample-name ${NORMAL} \
       --pass-only \
       --binding-threshold 500 \
       --top-score-metric median \
       --minimum-fold-change 1 \
       --normal-cov 5 \
       --tdna-cov 10 \
       --trna-cov 10 \
       --normal-vaf 0.02 \
       --tdna-vaf 0.25 \
       --trna-vaf 0.25 \
       --expn-val 1.0 \
       --maximum-transcript-support-level 1 \
       /vcf/${TUMOR}.vep.gx.cov.vcf \
       ${TUMOR} \
       ${HLA} \
       MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
       /pvac/
```
## Filter neoantigens
```
dir_MHC=${dir_PVAC}/${TUMOR}/MHC_Class_I
09_FilterEpitopes_new.sh ${dir_MHC} ${TUMOR}
dir_MHC=${dir_PVAC}/${TUMOR}/MHC_Class_II
09_FilterEpitopes_new.sh ${dir_MHC} ${TUMOR}
```
```
docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:3.0.1 \
        pvacseq binding_filter \
        --binding-threshold 99999999 \
        --minimum-fold-change 0 \
        --top-score-metric median \
        /mhc/${SAMPLE}.all_epitopes.tsv \
        /mhc/${SAMPLE}.all_epitopes.binding.tsv

docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:3.0.1 \
        pvacseq coverage_filter \
        --normal-cov 0 \
        --tdna-cov 0 \
        --trna-cov 0 \
        --normal-vaf 1 \
        --tdna-vaf 0 \
        --trna-vaf 0 \
        --expn-val 0 \
        /mhc/${SAMPLE}.all_epitopes.binding.tsv \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.tsv

docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:3.0.1 \
        pvacseq transcript_support_level_filter \
        --maximum-transcript-support-level 1 \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.tsv \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.transcript.tsv

docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:3.0.1 \
        pvacseq top_score_filter \
        --top-score-metric median \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.transcript.tsv \
        /mhc/${SAMPLE}.man_filter.tsv
```
