#!/bin/bash

# Find installation directory
install_dir=$(echo "$PATH" | tr : '\n' | grep focalscan)
if [ -z "$install_dir" ]
then
#    echo "The installation directory was not found in the PATH variable. Checking if the currect directory contains the necessary files..."
    if [ -f "./FocalScan.m" ];
    then
#        echo "Main file found. Temporarily adding current directory to PATH."
        export PATH="$PATH":"$(pwd)"
    else
        echo "The necessary files were not found in the current directory. Please manually add the directory countaining the FocalScan files to PATH by typing: export PATH=\$PATH:path_to/focalscan"
        exit $?
    fi

    install_dir=$(echo "$PATH" | tr : '\n' | grep focalscan)
    if [ -z "$install_dir" ]
    then
        echo "Failed adding current directory to PATH, for unknown reasons."
        exit $?
    fi
fi

current_dir=$(pwd)

argstring=$("$install_dir"/build_argstring.sh "$@")
#argstring=$(build_argstring.sh $@)

if [ -z "$argstring" ]
then
    echo
    echo "$(tput bold)USAGE:$(tput sgr0)
    standalone_peakdetection.sh report_file_path annot_file_path peak_level scorefield out_file

        $(tput bold)report_file_path:$(tput sgr0) path to main report file created by running FocalScan (ie. report.txt)
        $(tput bold)annot_file_path:$(tput sgr0) path to annotation file (ie. annotation/gencode17.bed, gene/tile ids must match those in the report file)
        $(tput bold)peak_level:$(tput sgr0) level of granularity at which to find peaks (0.1-1, where 1 is least granular)
        $(tput bold)scorefield:$(tput sgr0) metric in the report file to be used as basis for peak detection;
            (Valid options:
                $(tput smul)fs_hp:$(tput rmul) the standard FocalScan score (with focality filter)
                $(tput smul)fs:$(tput rmul) FocalScan score without focality filter
                $(tput smul)sum_cna_hp:$(tput rmul) summed copy number amplitudes, with focality filter 
                $(tput smul)sum_cna:$(tput rmul) summed copy number amplitudes, without focality filter 
                $(tput smul)spearman_corr:$(tput rmul) spearman correlation coefficient)
        $(tput bold)out_file:$(tput sgr0) name of output file (ie. peaks.txt)"
    echo
    echo "$(tput bold)EXAMPLE:$(tput sgr0)
    standalone_peakdetection.sh report.txt annotation/gencode17.bed 0.7 fs_hp ./peaks.txt"
    echo
else
    matlab -nodesktop -nosplash -r "addpath(genpath('$install_dir'));cd('$current_dir');try; standalone_peakdetection($argstring); catch ME; disp(getReport(ME,'extended','hyperlinks','off')); exit; end; exit"
fi
