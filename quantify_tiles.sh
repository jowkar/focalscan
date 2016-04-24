#!/bin/bash

# Assumes that the path to coverageBed is in the environment variable $PATH
# If not, download bedtools from https://github.com/arq5x/bedtools2/releases and add the following to ~/.bash_profile (or ~/.bash_rc, depending on system):
# export PATH=$PATH:path_to_bedtools/bin
# the load the new configuration:
# source ~/.bash_profile


bam_name=`echo "$1" | awk -F".bam" '{print $1}'`

if [ -z "$2" ]
then 
    annot=$2
else
    annot=./annotation/hg18_hg19_1kb_tiles.bed
fi

echo begin processing $1
# Quantify tiles
(coverageBed -split -counts -abam $1 -b $annot | cut -f4,5 > $bam_name.tile_counts) >& $bam_name.bedtools_log
echo $1: tiles quantified
