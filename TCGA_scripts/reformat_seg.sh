#!/bin/bash
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $1 | grep "TCGA-[A-Z0-9]\{2\}-[A-Z0-9]\{4\}-01" > grep.seg
cut -f2-6 grep.seg > SEG.txt
awk -F "-" '{print $1"-"$2"-"$3}' grep.seg > TC.txt
paste TC.txt SEG.txt > $1.reformatted.seg
rm TC.txt
rm grep.seg
rm SEG.txt
str=`head -n 1 $1`
sed "1s/^/$str\n/" $1.reformatted.seg > temp
mv temp $1.reformatted.seg
