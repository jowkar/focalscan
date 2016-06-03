#!/bin/bash

# Check if MATLAB is in PATH
matlab_path=$(which matlab)
if [ -z "$matlab_path" ]
then
    echo "Could not find the MATLAB executable. Please specify the path to MATLAB by typing:
    export PATH=\$PATH:path_to/matlab_directory/ 
    (for instance, export PATH=\$PATH:/Applications/MATLAB_R2016a.app/bin)"
    exit $?
fi

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
        echo "The necessary files were not found in the current directory. Please manually add the directory containing the FocalScan files to PATH by typing: export PATH=\$PATH:path_to/focalscan_directory"
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

if [ -z "$argstring" ]
then
    echo
    echo "FOCALSCAN: Rank genes based on coordinated focal changes in copy number and expression."
    echo
    echo "$(tput bold)USAGE:$(tput sgr0)
    focalscan.sh parameter_name1 parameter_value1 ... parameter_nameN parameter_valueN"

    echo
    echo "$(tput bold)EXAMPLE 1 (expression data in CSV format):$(tput sgr0)
    focalscan.sh expr_csv example_data/BRCA_expr.csv seg_file example_data/BRCA.seg annot_file annotation/gencode17_symbols.bed reportdir results"
    echo
    echo "$(tput bold)EXAMPLE 2 (separate read count files):$(tput sgr0)
    focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .gene_counts seg_file example_data/BRCA.seg annot_file annotation/gencode17.bed reportdir results"
echo
echo "$(tput bold)GENERAL USAGE INFO:$(tput sgr0)
    $(tput bold)-$(tput sgr0) The minimum required input is the path to expression data, copy number data and an annotation file.

    $(tput bold)-$(tput sgr0) The expression data can be given either as a CSV file (for gene-level analysis) or as a directory of separate read count files (gene-level and tile-level analysis). When the second option is used, the file extension of these files must also be given, and additionally an index file. The index file should contain two columns, where the first one includes the file names of the read count files (excluding the extension) and the second one includes the sample names used in the copy number data.

    $(tput bold)-$(tput sgr0) The copy number data needs to be provided in segmented format, in a file following the SEG format (https://www.broadinstitute.org/igv/SEG)

    $(tput bold)-$(tput sgr0) To perform a tile-level analysis, use the included file hg18_hg19_1kb_tiles.bed as annotation

    $(tput bold)-$(tput sgr0) For tile-level analysis, also remember to specify the parameter "optional_gene_annot" and give the path to a standard gene annotation file (BED format) in order to report which genes overlap each tile in the final peak report."
echo
echo "$(tput bold)PARAMETERS:$(tput sgr0)"
echo "
    [$(tput bold)Name:$(tput sgr0)]                    [$(tput smul)Default:$(tput rmul)]

    [Input]
    $(tput bold)expr_csv$(tput sgr0)*                    $(tput smul)''$(tput rmul)
    $(tput bold)seg_file$(tput sgr0)                     $(tput smul)''$(tput rmul)
    $(tput bold)annot_file$(tput sgr0)                   $(tput smul)''$(tput rmul)
    $(tput bold)expr_path$(tput sgr0)*                   $(tput smul)''$(tput rmul)
    $(tput bold)index_file$(tput sgr0)*                  $(tput smul)''$(tput rmul)
    $(tput bold)file_extension$(tput sgr0)*              $(tput smul)''$(tput rmul)
    $(tput bold)optional_gene_annot$(tput sgr0)          $(tput smul)''$(tput rmul)
    $(tput bold)fast_read$(tput sgr0)                    $(tput smul)0$(tput rmul)

    *Either 'expr_csv' or the combination of 'expr_path', 'index_file' and 'file_extension' should be used for expression data input.

    [Normalization]
    $(tput bold)normalization$(tput sgr0)                $(tput smul)percentile$(tput rmul) {$(tput smul)percentile$(tput rmul),$(tput smul)library_size$(tput rmul),$(tput smul)none$(tput rmul)}
    $(tput bold)percentile$(tput sgr0)                   $(tput smul)95$(tput rmul)

    [Score calculation]
    $(tput bold)window_size$(tput sgr0)                  $(tput smul)10e6$(tput rmul)
    $(tput bold)neutral_thresh$(tput sgr0)               $(tput smul)0.1$(tput rmul)
    $(tput bold)min_neutral$(tput sgr0)                  $(tput smul)20$(tput rmul)
    $(tput bold)pseudo_expr$(tput sgr0)                  $(tput smul)''$(tput rmul)
    $(tput bold)pseudo_expr_relative$(tput sgr0)         $(tput smul)10$(tput rmul)
    $(tput bold)max_nan$(tput sgr0)                      $(tput smul)0.1$(tput rmul)

    [Output and peak detection]
    $(tput bold)peak_level$(tput sgr0)                   $(tput smul)0.6$(tput rmul) {$(tput smul)0.0$(tput rmul)-$(tput smul)1.0$(tput rmul)}
    $(tput bold)reportdir$(tput sgr0)                    $(tput smul)'.'$(tput rmul)
    $(tput bold)only_focal$(tput sgr0)                   $(tput smul)''$(tput rmul)
    $(tput bold)scorefield$(tput sgr0)                   $(tput smul)fs_hp$(tput rmul)

$(tput bold)PARAMETER DESCRIPTIONS:$(tput sgr0)
    $(tput bold)annot_file:$(tput sgr0) Gene annotation or tile definition file in .bed format.
    $(tput bold)expr_csv:$(tput sgr0) Path to a CSV file containing unnormalized expression data. Columns are expected to correspond to samples and rows to genes. The columns should be titled with sample IDs. (Only for gene-level analysis)
    $(tput bold)expr_path:$(tput sgr0) Path to directory containing files with gene or tile level count data for all samples (given in separate files)
    $(tput bold)fast_read:$(tput sgr0) When set to 1 and separate read count files are used, will assume that all files have identical first columns (gene IDs) in order to speed up reading of these.
    $(tput bold)file_extension:$(tput sgr0) The file extension of the gene or tile-level expression files.
        # {expr path, index file, file extension}: Need to be specified together.
    $(tput bold)index_file:$(tput sgr0) File that links expression data files to sample IDs
    $(tput bold)max_nan:$(tput sgr0) Maximum proportion of missing values to accept for a given gene/tile
    $(tput bold)neutral_thresh:$(tput sgr0) Absolute copy number amplitude threshold for defining
    $(tput bold)neutral samples$(tput sgr0)
    $(tput bold)normalization:$(tput sgr0) The normalization mode to employ
    $(tput bold)only_focal:$(tput sgr0) When set to 1, will avoid additional calculation of scores without the focality filter (will speed up execution).
    $(tput bold)optional_gene_annot:$(tput sgr0) When tile-level analysis is performed, providing a gene annotation via this option will enable annotation of the reported peak tiles with respect to overlapping genes.
    $(tput bold)peak_level:$(tput sgr0) Sets the granularity of the peak detection method. A high value will cause only the most prominent peaks to be reported. A low value will cause additional, less prominent, peaks to be reported.
    $(tput bold)percentile:$(tput sgr0) The percentile to use when percentile normalization is employed. For instance, '95' will normalize to the median of the top 5 percent most highly expressed genes in each sample
    $(tput bold)pseudo_expr:$(tput sgr0) Pseudo expression value to add (needed to avoid division with zero when calculating ratios)
    $(tput bold)pseudo_expr_relative:$(tput sgr0) The pseudo expression value can be specified relative to the median of all non-zero expression values. This parameter defined the relation between the pseudo count and this median. For instance, a value of 10 sets the pseudo count to 10 times the median.
    $(tput bold)reportdir:$(tput sgr0) Directory in which to store output files
    $(tput bold)scorefield:$(tput sgr0) The metric to use as basis for peak detection.
    $(tput bold)seg_file:$(tput sgr0) File containing segmented copy number data for all samples 
    $(tput bold)window_size:$(tput sgr0) Window size used by the focality filter 

"
else
    matlab -nodesktop -nosplash -r "addpath(genpath('$install_dir'));cd('$current_dir');try; FocalScan($argstring); catch ME; disp(getReport(ME,'extended','hyperlinks','off')); exit; end; exit"
fi
