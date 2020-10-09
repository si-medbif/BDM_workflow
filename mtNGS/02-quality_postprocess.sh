#!/bin/bash

#Collect alignment summary metrics
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do gatk CollectAlignmentSummaryMetrics -R /penguin/resources/rCRS/rCRS.fasta -I /turtle/BSL/mon/aligned/$PREFIX.readgroup.bam -O /turtle/BSL/mon/stats/$PREFIX.alignment.metrics; done

#Collect the total number of bases
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do zcat /turtle/BSL/mon/fastq/${PREFIX}_L001_R*_001.fastq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c > stat.tmp && awk '{print "TotaL Number of Bases = " $1}' stat.tmp >> /turtle/BSL/mon/stats/$PREFIX.summary.metrics && rm stat.tmp; done

#Collect the total number of fragments
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do grep -w 'PAIR' /turtle/BSL/mon/stats/$PREFIX.alignment.metrics | awk '{print "Total Number of Fragments = " $2}' > /turtle/BSL/mon/stats/$PREFIX.summary.metrics; done

#Collect the total number of high quality fragments
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do awk '(NR > 389 && NR < 399) {print $2}' /turtle/BSL/mon/stats/${PREFIX}_L001_R1_001_fastqc/fastqc_data.txt > stat.tmp && awk '(NR > 389 && NR < 399) {print $2}' /turtle/BSL/mon/stats/${PREFIX}_L001_R2_001_fastqc/fastqc_data.txt >> stat.tmp && awk '{total += $1} END {print "Total Number of High Quality Reads= " total}' stat.tmp >> /turtle/BSL/mon/stats/$PREFIX.summary.metrics && rm stat.tmp; done

#Calculate average fragment length
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do zcat /turtle/BSL/mon/fastq/${PREFIX}_L001_R*_001.fastq.gz | awk '{if(NR%4==2) {count++; bases += length}} END {print "Average Fragement Length = " bases/count}' >> /turtle/BSL/mon/stats/$PREFIX.summary.metrics; done

#Calculate average depth across mitochondrial genome and mitochondrial coverage
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do samtools depth /turtle/BSL/mon/aligned/$PREFIX.readgroup.bam | awk '{total+=$3} END {print "Average Depth = " total/16569 "\nMitochondrial Coverage = " NR/16569*100 "%"}' >> /turtle/BSL/mon/stats/$PREFIX.summary.metrics; done

#Calculate Mapping Quality scores
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do bam stats --in /turtle/BSL/mon/aligned/$PREFIX.readgroup.bam --cBaseQC /turtle/BSL/mon/stats/$PREFIX.tmp && awk '{total += $14} END {print "Mapping Quality = " total/NR}' /turtle/BSL/mon/stats/$PREFIX.tmp >> /turtle/BSL/mon/stats/$PREFIX.summary.metrics && rm /turtle/BSL/mon/stats/$PREFIX.tmp; done
