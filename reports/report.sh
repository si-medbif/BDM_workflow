#! /bin/bash

#set -euo pipefail

IFS=$'\n'

SOFTWARE=$1
SAMPLES=$2

for input in $(cat ${SAMPLES}) ;
do
	IFS=' '
	read -a strarr <<< "$input"
	./${SOFTWARE} ${strarr[0]} ${strarr[1]}
	IFS=$'\n'
done
