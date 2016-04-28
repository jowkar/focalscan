#!/bin/bash

# Find installation directory
install_dir=`echo $PATH | tr : '\n' | grep focalscan`
if [ -z "$install_dir" ]
then
    echo "The installation directory was not found in the PATH variable. Checking if the currect directory contains the necessary files"
    main_file=`ls | grep FocalScan.m`
    echo $main_file
    if [ -z "$main_file" ]
    then
        echo "The necessary files were not found in the current directory. Please manually add the directory countaining the FocalScan files to PATH by typing: export PATH=\$PATH:path_to/focalscan"
        exit $?
    else
        echo "Main file found. Adding current directory to PATH"
        export PATH=$PATH:`pwd`
    fi
    install_dir=`echo $PATH | tr : '\n' | grep focalscan`
    if [ -z "$install_dir" ]
    then
        echo "Failed adding current directory to PATH, for unknown reasons"
        exit $?
    fi
fi

current_dir=`pwd`

argstring=`$install_dir/build_argstring.sh $@`
#argstring=`build_argstring.sh $@`

matlab -nodesktop -nosplash -r "addpath(genpath('$install_dir'));cd('$current_dir');try; standalone_peakdetection($argstring); catch ME; disp(getReport(ME,'extended','hyperlinks','off')); exit; end; exit"
