#!/bin/bash

#Align reads to a reference genome
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do bwa mem -t 64 /penguin/resources/rCRS/rCRS.fasta /turtle/BSL/mon/fastq/${PREFIX}_L001_R1_001.fastq.gz /turtle/BSL/mon/fastq/${PREFIX}_L001_R2_001.fastq.gz | samtools view -Shb - | samtools sort - -o /turtle/BSL/mon/aligned/$PREFIX.sorted.bam; done

#Add read groups
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do gatk AddOrReplaceReadGroups -I /turtle/BSL/mon/aligned/$PREFIX.sorted.bam -O /turtle/BSL/mon/aligned/$PREFIX.readgroup.bam -ID $PREFIX -LB Mitochondria -PL illumina -PU NA -SM Mon; done

#Index BAM files
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do samtools index /turtle/BSL/mon/aligned/$PREFIX.readgroup.bam; done
