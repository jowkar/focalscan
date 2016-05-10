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
#argstring="$argstring,'current_dir','$current_dir'"

#matlab -nodesktop -nosplash -r "FocalScan($argstring);exit"

if [ -z "$argstring" ]
then
    echo
    echo "FOCALSCAN: Rank genes based on coordinated focal changes in copy number and expression."
    echo
    echo "$(tput bold)USAGE:$(tput sgr0)
    focalscan.sh parameter_name1 parameter_value1 ... parameter_nameN parameter_valueN

    [$(tput bold)Parameter:$(tput sgr0)]                    [$(tput smul)Default:$(tput rmul)]
    
    [Basic options]
	$(tput bold)expr_csv$(tput sgr0)                     $(tput smul)''$(tput rmul)
	$(tput bold)seg_file$(tput sgr0)                     $(tput smul)''$(tput rmul)
	$(tput bold)annot_file$(tput sgr0)                   $(tput smul)''$(tput rmul)
	
    [Alternative input options for expression data]
	$(tput bold)expr_path$(tput sgr0)                    $(tput smul)''$(tput rmul)
	$(tput bold)index_file$(tput sgr0)                   $(tput smul)''$(tput rmul)
	$(tput bold)file_extension$(tput sgr0)               $(tput smul)''$(tput rmul)
	$(tput bold)expr_ratio_csv$(tput sgr0)               $(tput smul)''$(tput rmul)
	
    [Additional options]
	$(tput bold)window_size$(tput sgr0)                  $(tput smul)10e6$(tput rmul)
	$(tput bold)neutral_thresh$(tput sgr0)               $(tput smul)0.1$(tput rmul)
	$(tput bold)min_neutral$(tput sgr0)                  $(tput smul)20$(tput rmul)
	$(tput bold)pseudo_expr$(tput sgr0)                  $(tput smul)''$(tput rmul)
	$(tput bold)pseudo_expr_relative$(tput sgr0)         $(tput smul)10$(tput rmul)
	$(tput bold)max_nan$(tput sgr0)                      $(tput smul)0.1$(tput rmul)
	$(tput bold)reportdir$(tput sgr0)                    $(tput smul)'.'$(tput rmul)
	$(tput bold)normalization$(tput sgr0)                $(tput smul)percentile$(tput rmul) {$(tput smul)percentile$(tput rmul),$(tput smul)library_size$(tput rmul),$(tput smul)none$(tput rmul)}
	$(tput bold)percentile$(tput sgr0)                   $(tput smul)95$(tput rmul)
	$(tput bold)optional_gene_annot$(tput sgr0)          $(tput smul)''$(tput rmul)
	$(tput bold)peak_level$(tput sgr0)                   $(tput smul)0.6$(tput rmul) {$(tput smul)0.0$(tput rmul)-$(tput smul)1.0$(tput rmul)}
	$(tput bold)only_focal$(tput sgr0)                   $(tput smul)''$(tput rmul)
	$(tput bold)scorefield$(tput sgr0)                   $(tput smul)fs_hp$(tput rmul)

"
    echo "$(tput bold)EXAMPLE:$(tput sgr0)
    focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .gene_counts seg_file example_data/BRCA.seg annot_file annotation/gencode17.bed reportdir results"
    echo
else
    matlab -nodesktop -nosplash -r "addpath(genpath('$install_dir'));cd('$current_dir');try; FocalScan($argstring); catch ME; disp(getReport(ME,'extended','hyperlinks','off')); exit; end; exit"
fi
