#!/bin/bash

PROJECT=$1

mkdir -p ${PROJECT}

cp RNA/fastqc.sh ${PROJECT}/00_fastqc.sh
cp RNA/align.sh ${PROJECT}/01_align.sh
cp RNA/count.sh ${PROJECT}/02_count.sh
cp reports/count2tpm.py ${PROJECT}/03_count2tpm.py
