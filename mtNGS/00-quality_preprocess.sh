#!/bin/bash

#Check sample quality
fastqc -o /turtle/BSL/mon/stats --noextract -t 64 /turtle/BSL/mon/fastq/*.fastq.gz
