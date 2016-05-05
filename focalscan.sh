#!/bin/bash
#argstring=""
#for var in "$@"
#do
#    argstring = $argstring,\'$var\'
#done


# Find installation directory
install_dir=$(echo "$PATH" | tr : '\n' | grep focalscan)
if [ -z "$install_dir" ]
then
    echo "The installation directory was not found in the PATH variable. Checking if the currect directory contains the necessary files..."
    if [ -f "./FocalScan.m" ];
    then
        echo "Main file found. Temporarily adding current directory to PATH"
        export PATH="$PATH":"$(pwd)"
    else
        echo "The necessary files were not found in the current directory. Please manually add the directory countaining the FocalScan files to PATH by typing: export PATH=\$PATH:path_to/focalscan"
        exit $?
    fi

    install_dir=$(echo "$PATH" | tr : '\n' | grep focalscan)
    if [ -z "$install_dir" ]
    then
        echo "Failed adding current directory to PATH, for unknown reasons"
        exit $?
    fi
fi

current_dir=$(pwd)

argstring=$("$install_dir"/build_argstring.sh "$@")
#argstring="$argstring,'current_dir','$current_dir'"

#matlab -nodesktop -nosplash -r "FocalScan($argstring);exit"

matlab -nodesktop -nosplash -r "addpath(genpath('$install_dir'));cd('$current_dir');try; FocalScan($argstring); catch ME; disp(getReport(ME,'extended','hyperlinks','off')); exit; end; exit"
