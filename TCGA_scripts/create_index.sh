#!/bin/bash
summary_file=$1
out_file=$2
run=1
if [ -z "$summary_file" ]
then
    echo "The first argument should be the path to a cghub summary file"
    run=0
fi
if [ -z "$out_file" ]
then
    echo "The second argument should be the output file name"
    run=0
fi
if [ "$run" -eq 1 ]
then
    cat $summary_file | awk '{print $2"\t"$14}' | tail -n +2 | awk -F".bam" '{print $1}' | awk '{print $2"\t"$1}' > $out_file
    cat $out_file | awk '{print $2}' | awk -F"-" '{print $1"-"$2"-"$3}' > tmp.txt
    awk '{print $1}' $out_file > tmp2.txt
    paste tmp2.txt tmp.txt > $out_file
    rm tmp.txt
    rm tmp2.txt
fi
