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
        echo "Failed adding current directory to PATH, for unknown reasons. Please manually add the directory countaining the FocalScan files to PATH by typing: export PATH=\$PATH:path_to/focalscan"
        exit $?
    fi
fi

current_dir=$(pwd)

argstring=$("$install_dir"/build_argstring.sh "$@")
#argstring=$(build_argstring.sh $@)

if [ -z "$argstring" ]
then
    echo
    echo "Add gene IDs to peak report file created by FocalScan. Useful if a tile-level analysis was performed without providing an optional gene annotation."
    echo
    echo "$(tput bold)USAGE:$(tput sgr0)
    annotate_peaks.sh peaks_file_path annot_file_path out_file

        $(tput bold)peak_file_path:$(tput sgr0) path to main report file created by running FocalScan (ie. peaks.txt)
        $(tput bold)annot_file_path:$(tput sgr0) path to annotation file (ie. annotation/gencode17.bed)
        $(tput bold)out_file:$(tput sgr0) name of output file (ie. peaks_annotated.txt)"
    echo
    echo "$(tput bold)EXAMPLE:$(tput sgr0)
    annotate_peaks.sh peaks.txt annotation/gencode17.bed ./peaks_annotated.txt"
    echo
else
    matlab -nodesktop -nosplash -r "addpath(genpath('$install_dir'));cd('$current_dir');try; annotate_peaks($argstring); catch ME; disp(getReport(ME,'extended','hyperlinks','off')); exit; end; exit"
fi
