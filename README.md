# About

FocalScan identifies genomic regions where many tumors show simultaneous increases in DNA copy-number amplitude (CNA) and RNA expression (or conversely for DNA deletions). Empirically, many important oncogenes show this pattern of alteration. The FocalScan score is based on the dot product of CNA and RNA changes. This puts equal weight to the two variables, but requires coordinated changes in both to achieve a positive score.

Examples (for a given genomic position):

Some tumors show both elevated CNA and RNA levels -> medium score
Many tumors show both elevated CNA and RNA levels -> high score
Many tumors show elevated CNA and highly increased RNA -> very high score Many tumors show elevated CNA, but RNA is unchanged -> neutral score

Regions with coordinated CNA and RNA reduction will also score favorably.

FocalScan computes two basic statistics: One is calculated as described above. The other is based on ’high-pass filtered’ CNA data, where large (>10 Mbp) segments, such as arm-level events, are effectively subtracted, leaving only the focal/small alterations. Apart from that, the same scoring method is used in both cases. The second method can be very sensitive at identifying genes of interest in focally altered regions.

Importantly, FocalScan can also be run in a “non gene-centric” fashion: The genome is scanned at high (500 nt) resolution by dividing chromosomes into small (1000 nt) overlapping tiles. As such, it does not care about preconceptions about gene locations. RNA-seq data is used to quantify transcription and scores are computed for each tile. This makes FocalScan suitable for identifying e.g. novel non-coding RNAs that are altered in tumors.


# System requirements

Linux/Mac/Windows, 8 GB RAM preferred for gene-based analysis, >30 GB RAM preferred for tile-based analysis. The included shell scripts (.sh) are not supported on Windows (only usage from within MATLAB possible).

Required software:

- MATLAB (unless the compiled executable is uesd)

Recommeded for pre-processing expression data:

- bedtools (http://bedtools.readthedocs.org/en/latest/)
- samtools (http://www.htslib.org)
- HTseq (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)

Recommended for visualizing the results:

- IGV (https://www.broadinstitute.org/igv/)


# 1. Installation

There are three main ways in which this tool can be run:

- Using the shell scripts (Linux/Mac, MATLAB 2015b or later)
- Using the compiled executable (Linux/Mac)
- From within the MATLAB environment (Linux/Mac/Windows, MATLAB 2015b or later)

Using the shell scripts:
------------------------

1. Open the terminal application and download the scripts:
```git clone git@github.com:jowkar/focalscan.git```

2. Setup environment variables (in the following, substitute path to focalscan for the directory where the files were downloaded (for instance, ∼/focalscan) and path to matlab for the directory where MATLAB is installed (for in- stance, ∼/bin/MATLAB 2015b/bin/)):

    ```shell
    export PATH=$PATH:path_to_focalscan
    export PATH=$PATH:path_to_matlab
    ```

    Example:

    ```shell
    export PATH=$PATH:~/focalscan
    export PATH=$PATH:~/bin/MATLAB_R2016a/bin
    ```

    The above paths will differ depending on system and where the respective applications were installed by the user.

3. *Optional:* To avoid having to type the commands in the previous step every time the tool is run, add those commands to the ∼/.bashrc or ∼/.bash profile file (depending on the operating system).

4. Add permissions to execute the shell scripts:

    ```shell
    chmod +x *.sh
    ```

5. Download and unzip example data to test the installation with, available at the following links:

    _Annotation files:_

    ```shell
    https://drive.google.com/open?id=0B_52viSz8FLNeUlPZ0c3akU1ZlE
    ```

    _Test data:_

    ```shell
    https://drive.google.com/open?id=0B_52viSz8FLNWFRvaUNyMGhCbHM
    ```

6. Test the installation (gene-level analysis on breast cancer data from TCGA) (all in one line):

    ```shell
    focalscan.sh expr_csv ./example_data/BRCA_expr.csv seg_file ./example_data/BRCA_cna.seg annot_file ./annotation/gencode17_symbols.bed reportdir test_gene
    ```

7. Inspect the output:

    ```shell
    cat test_gene/peaks.txt
    ```

    This should list the top ranking genes:

    ```
    Id  Score   Sum_CNA_HP  Chr Start   Stop
    ERBB2   942.124267578125    262.793273925781    37844167    37886679
    CCND1   447.170532226562    220.362243652344    69455855    69469242
    WHSC1L1 440.08837890625 202.27294921875 chr8    38132544    38239790
    EGFR    129.50341796875 21.137996673584 chr7 55086714   55324313
    TRAF4   117.803070068359    80.5899887084961    chr17   27071002    27077974
    IGF1R   116.806274414062    26.7439022064209    chr15   99192200    99507759
    FGFR2   96.7242736816406    24.9825992584229    chr10   123237848   123357972
    CCNE1   83.940673828125 27.0290393829346    chr19   30302805    30315215
    PHGDH   80.3051528930664    20.5632972717285    chr1    120202421   120286838
    ```

    Also available in the in the output are a full report with detailed statistics and .wig files that can be visualized with IGV, containing full tracks with scores, copy number amplitudes and mean expression levels of all genes.

Mac/Linux, compiled executable
------------------------------

1. Download the file titled FocalScan compiled Linux.zip or FocalScan compiled Mac.zip.

2. Open the terminal application and unarchive the file:

    ```shell
    unzip FocalScan_compiled_Linux.zip
    ```

3. Install the MATLAB runtime (it is important that the runtime is the correct version, in this case v901):

    ```shell
    cd FocalScan_compiled_Linux
    ./Installer.install
    ```
    or

    ```shell
    cd FocalScan_compiled_Mac
    ./Installer.app
    ```
    This should bring up a window for downloading and installing the runtime. Make a note of where FocalScan is installed and where the runtime is installed (such as /Applications/FocalScan/application and ∼/bin/MCR/v901, re- spectively, or any other directories chosen). The path to the runtime directory (from now on referred to as “MCR root”) and the path to FocalScan (from now on referred to as “path to focalscan”) will have to be specified later when running the program. Note also that the program files will be located in a subdirectory named “application”. I should also be noted that remote installation on a Linux server was observed to occassionally fail if X forwarding was not used (a bug in the MATLAB installer program).

    4. Setup environment variables (in the following, substitute path to focalscan for the directory where the files were downloaded (for instance, ∼/focalscan):

    ```shell
    export PATH=$PATH:path_to_focalscan
    ```

    Example:

    ```shell
    export PATH=$PATH:/Applications/FocalScan/application
    ```

5. Download and unzip example data to test the installation with, available at the following links:

    _Annotation files:_

    ```shell
    https://drive.google.com/open?id=0B_52viSz8FLNeUlPZ0c3akU1ZlE
    ```

    _Test data:_

    ```shell
    https://drive.google.com/open?id=0B_52viSz8FLNWFRvaUNyMGhCbHM
    ```

6. Test the installation (gene-level analysis on breast cancer data from TCGA). In the following command, substitute “MCR root” for the installation direc- tory of the MATLAB compiler runtime and path to focalscan for the installation directory of FocalScan (obtained from step 3 above) (all in one line):

    ```shell
    focalscan_compiled.sh MCR_root expr_csv ./example_data/BRCA_expr.csv seg_file ./example_data/BRCA_cna.seg annot_file ./annotation/gencode17_symbols.bed reportdir test_gene
    ```

    This might take up to 40 minutes or slightly longer, depending on processor speed and available memory.

7. Inspect the output:

    ```shell
    cat test_gene/peaks.txt
    ```

Any platform, usage from within the MATLAB environment
------------------------------------------------------

1. Open the terminal application and download the scripts (assuming that git is installed, otherwise just download the zip file):

    ```shell
    git clone git@github.com:jowkar/focalscan.git
    ```

2. Download and unzip example data to test the installation with, available at the following links:

    _Annotation files:_

    ```
    https://drive.google.com/open?id=0B_52viSz8FLNeUlPZ0c3akU1ZlE
    ```

    _Test data:_

    ```
    https://drive.google.com/open?id=0B_52viSz8FLNWFRvaUNyMGhCbHM
    ```

3. Open MATLAB and add the scripts to the MATLAB path (will have to re- peated each time MATLAB is opened, unless the startup settings are also changed):

    ```MATLAB
    addpath(genpath('path_to_focalscan'))
    ```

    Example: 

    ```MATLAB
    addpath(genpath('~/focalscan'))
    ```

    (If focalscan was downloaded to the home directory on Linux/Mac)

4. Test the installation (gene-level analysis on breast cancer data from TCGA):

    ```MATLAB
    FocalScan.sh('expr_csv','example_data/BRCA_expr.csv','seg_file','example_data/BRCA_cna.seg','annot_file','annotation/gencode17_symbols.bed','reportdir','test_gene')
    ```

    This might take up to 40 minutes or slightly longer, depending on proces- sor speed and available memory. The results will be saved to the directory “test gene”.

# 2. Input files

RNA-seq data should be provided as read counts per gene/tile. Gene-level read counts can be obtained with HTSeq (or equivalent) and tile-level counts with coverageBed (for instance using the included script "quantify_tiles.sh").

Two FocalScan expression data input options exist:
--------------------------------------------------

1. A CSV-file with samples as columns and genes as rows.

2. A combination of the following:
    - The path to a directory containing separate read count files for each sample
    - An index file mapping the file names in this directory to sample names
    - The file extension (for instance ".gene_counts") of the read count files in this directory

    The index file is a tab-delimited file connecting read count files (excluding the file extension) and tumor ID's:
    ```
    RNA-seqfilename1    tumorID1
    RNA-seqfilename2    tumorID2
    ```

Segmented copy-number data should be provided as a single .seg file with data for all tumors.
---------------------------------------------------------------------------------------------

NOTE:
-----
- Tumor ID's in the index file need to match those in .seg file.

- Ensure that both the expression data files and the .seg file refers to the same genome assembly (e.g. hg18 and hg19).

- Only tumors - do not include normal samples.

- Likely, germline variants need to be filtered out from the copy number data beforehand.


# 3. Running FocalScan


## 1) Pre-process the data:

Gene-level analysis:
--------------------

Use HTSeq or equivalent to obtain gene read counts for each sample.

Tile-level analysis:
--------------------
```shell
./quantify_tiles.sh <BAM_file> annotation/hg18_hg19_1kb_tiles.bed
```
This needs to be done for each .bam file, preferably in parallell to speed up processing.
A file with RNA-seq read counts (\*.tile_counts) for each genomic tile will be generated.

NOTE: Your .bam files may use chromosome names formatted as e.g. 'chr1' or simply '1'. In
the latter case, instead use the hg18_hg19_1kb_tiles_nochr.bed file.

NOTE: You may consider pre-filtering your .bam files to only consider uniquely mapped/high
quality reads (e.g. quality 255 only for TopHat alignments, by running 'samtools view -b
-q255 in.bam > out.bam').

## 2) Run FocalScan:
-----------------


- Using the shell scripts:
```shell
./focalscan.sh <parameter1> <parameter1_value> ... <parameterN> <parameterN_value>
```
Example:

Using CSV input:
```shell
./focalscan.sh expr_csv example_data/BRCA_expr.csv seg_file example_data/BRCA_cna.seg annot_file annotation/gencode17_symbols.bed
```

Using separate read count files for each sample:
```shell
./focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .gene_counts annot_file annotation/gencode17.bed seg_file example_data/BRCA_cna.seg
```
Tile-level analysis:
```shell
./focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .tile_counts annot_file annotation/gencode17.bed seg_file example_data/BRCA_cna.seg
```
- From within MATLAB:
```MATLAB
FocalScan('parameter1','parameter1_value',...,'parameterN','parameterN_value')
```
Example:
```MATLAB
FocalScan('expr_csv','example_data/BRCA_expr.csv','seg_file','example_data/BRCA_cna.seg','annot_file','annotation/gencode17_symbols.bed')
```

- With the compiled executable (no MATLAB installation)

```shell
focalscan_compiled.sh <MCR_root> <parameter1> <parameter1_value> ... <parameterN> <parameterN_value>
```

where <MCR_root> is the path to the MATLAB runtime (v901, downloaded by the included installer).

See the manual for more information about available input options

Several output files will be generated (use the 'reportdir' parameter to specify where to write these files. '.' is default):

- report.txt: full report with statistics for each gene/tile
- peaks.txt: top ranked genes/tiles, based on peak detection
- log.txt: log file
- rna.wig: mean expression
- score_hp.wig: standard FocalScan score (with focality filter)
- score.wig: FocalScan score without focality filter
- sum_cna_hp.wig: Summed copy number amplitudes, with focality filter
- sum_cna.wig: Sumemd copy number amplitudes, without focality filter

The WIG files can be visualized with IGV

Standalone peak detection:
--------------------------

After running FocalScan, peak detection may be re-run if desired.

Either use the shell script:

```standalone_peakdetection.sh report_file_path annot_file_path peak_level scorefield out_file```

Valid options for the "scorefield" parameter (the metric to use as basis for peak detection) are:

- "fs_hp": the standard FocalScan score (with focality filter)
- "fs": FocalScan score without focality filter
- "sum_cna_hp": summed copy number amplitudes, with focality filter
- "sum_cna": summed copy number amplitudes, without focality filter
- "spearman_corr": spearman correlation coefficient

(Assuming that all of the above scores are present in the report.txt file, as is the case by default.)

Example:

```shell
standalone_peakdetection.sh example_data/test_CSV/report.txt annotation/gencode17_symbols.bed 0.7 fs_hp ./new_peaks.txt
```

Or the MATLAB function:

```shell
standalone_peakdetection(report_file_path,annot_file_path,peak_level,scorefield,out_file)
```

Or the compiled version 

```shell
standalone_peakdetection_compiled.sh <MCR_root> example_data/test_CSV/report.txt annotation/gencode17_symbols.bed 0.7 fs_hp ./new_peaks.txt
```

Annotating a tile peak report:
------------------------------

By default, the program will not write the IDs of genes overlapping the tiles in the output report unless a separate gene annotation file was added with the parameter "optional_gene_annot". To add this information afterwards the script "annotate_peaks.sh" can be used:

```shell
annotate_peaks.sh peak_file_path annot_file_path out_file
```

Example:

```shell
annotate_peaks.sh peaks.txt annotation/gencode17_symbols.bed peaks_annotated.txt
```

The compiled version is named "annotate_peaks_compiled.sh".
