# About

FocalScan identifies genomic regions where many tumors show simultaneous increases in DNA copy-number amplitude (CNA) and RNA expression (or conversely for DNA deletions). Empirically, many important oncogenes show this pattern of alteration. The FocalScan score is based on the dot product of CNA and RNA changes. This puts equal weight to the two variables, but requires coordinated changes in both to achieve a positive score.

Examples (for a given genomic position):

- Some tumors show both elevated CNA and RNA levels -> medium score
- Many tumors show both elevated CNA and RNA levels -> high score
- Many tumors show elevated CNA and highly increased RNA -> very high score 
- Many tumors show elevated CNA, but RNA is unchanged -> neutral score

Regions with coordinated CNA and RNA reduction will also score favorably.

FocalScan computes two basic statistics: One is calculated as described above. The other is based on ’high-pass filtered’ CNA data, where large (>10 Mbp) segments, such as arm-level events, are effectively subtracted, leaving only the focal/small alterations. Apart from that, the same scoring method is used in both cases. The second method can be very sensitive at identifying genes of interest in focally altered regions.

Importantly, FocalScan can also be run in a “non gene-centric” fashion: The genome is scanned at high (500 nt) resolution by dividing chromosomes into small (1000 nt) overlapping tiles. As such, it does not care about preconceptions about gene locations. RNA-seq data is used to quantify transcription and scores are computed for each tile. This makes FocalScan suitable for identifying e.g. novel non-coding RNAs that are altered in tumors.


# System requirements

Linux/Mac/Windows, 8 GB RAM preferred for gene-based analysis, >30 GB RAM preferred for tile-based analysis. The included shell scripts (.sh) are not supported on Windows (where only usage from within MATLAB possible).

**Note:** This is a command line based application.

**Operating system:**
- _Standard shell scripts, with MATLAB installed_ **(preferred)**:
    - Any operating system that supports MATLAB R2013b or newer. (Tested on OSX 10.11, CentOS 6.6, 6.7, Fedora 20, Windows 7)
- _Compiled executable_:
    - Mac: **OSX 10.11** (El Capitan). 
        - Not guarateed to work on older versions (in theory 10.10 should work as well, although this has not been tested). 10.8 (Mountain Lion) does **not** work (the v901 MATLAB runtime is incompatible).
    - Linux: CentOS 6.6, 6.7 tested
        - Ubuntu 14.04 and Fedora 20 have also been tested successfully. However, the following issues were encountered (which are operating system and/or MATLAB specific bugs unrelated to FocalScan itself):
            - Ubuntu 14.04: the installer appeared to get stuck while downloading the runtime from MathWorks, but clicking "Finish" after a while completed the installation sucessfully despite this.
            - Fedora 20: may have issues connecting to the MathWorks servers

**Hardware:**
- Gene-level analysis:
    - 8 GB RAM or more, depending on the number of samples analyzed
- Tile-level analysis:
    - 30 GB RAM or more, depending on the number of samples

**Required software:**

- MATLAB R2013b to R2016a (unless the compiled executable is used)

**Recommeded for pre-processing expression data:**

- bedtools 2.21.0 (http://bedtools.readthedocs.org/en/latest/)
- samtools 1.1 (http://www.htslib.org)
- HTseq 0.6.1 (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)

**Recommended for visualizing the results:**

- IGV 2.3.60 (https://www.broadinstitute.org/igv/)

# Additional files

- Annotation: https://drive.google.com/file/d/0B_52viSz8FLNTk4xb2VuY3lSQk0/view?usp=sharing
- Example data: https://drive.google.com/file/d/0B_52viSz8FLNX2lZVF9SOTFMQzg/view?usp=sharing
- Compiled Linux version: https://drive.google.com/file/d/0B_52viSz8FLNVmk2RUJOeU82M2M/view?usp=sharing
- Compiled Mac version: https://drive.google.com/file/d/0B_52viSz8FLNc2VvTERsUWVyVUU/view?usp=sharing

# 1. Installation

There are three main ways in which this tool can be run:

- **If MATLAB is installed:** 
    - Use the shell scripts ("focalscan.sh") available here at GitHub (Linux/Mac, MATLAB 2013b or later) **(preferred)**
    - From within the MATLAB environment (the function "FocalScan") (Linux/Mac/Windows, MATLAB 2013b or later)
- **If MATLAB is not installed:** 
    - Use the compiled executable (the script "focalscan", without ".sh", provided by the installer) (Linux/Mac)
        - See also the section on system requirements

**Note:** The manual contains a section detailing common errors and warnings that may be referred to if anything goes wrong.

**Note:** Using the standard shell scripts (ie. focalscan.sh) and a MATLAB installation is preferred, since the compiled executables are somewhat sensitive to differences between operating system versions, due to how the MATLAB runtime and its installer operates. For instance, Mac OSX 10.8 cannot execute the v901 runtime and some Linux versions have issues connecting to the MathWorks servers.

Using the shell scripts ("focalsan.sh"):
------------------------

1. Open the terminal application and download the scripts:
```git clone https://github.com/jowkar/focalscan.git```
(or just download the ZIP archive from https://github.com/jowkar/focalscan/archive/master.zip)

2. Setup environment variables (in the following, substitute path to focalscan for the directory where the files were downloaded (for instance, ∼/focalscan) and path to matlab for the directory where MATLAB is installed (for instance, ∼/bin/MATLAB 2013b/bin/)):

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

3. **Optional:** To avoid having to type the commands in the previous step every time the tool is run, add those commands to the ∼/.bashrc or ∼/.bash profile file (depending on the operating system).

4. Add permissions to execute the shell scripts:

    ```shell
    chmod +x *.sh
    ```

5. Download and unzip example data to test the installation with, available at the following links:

    **Annotation files:**

    - https://drive.google.com/file/d/0B_52viSz8FLNTk4xb2VuY3lSQk0/view?usp=sharing

    **Test data:**

    - https://drive.google.com/file/d/0B_52viSz8FLNX2lZVF9SOTFMQzg/view?usp=sharing

6. Test the installation (gene-level analysis on breast cancer data from TCGA) (all in one line):

    ```shell
    focalscan.sh expr_csv example_data/BRCA_expr.csv seg_file example_data/BRCA.seg annot_file annotation/gencode17_symbols.bed reportdir reports
    ```

    If the following occurs, double check that the correct input file paths were given:
    ```shell
    Invalid file identifier.  Use fopen to generate a valid file identifier.
    ```

7. Inspect the output:

    ```shell
    cat reports/peaks.txt
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

Mac/Linux, compiled executable (started by the script "focalscan", without ".sh")
------------------------------

1. Download the file Installer_Linux.install.zip or Installer_Mac.app.zip:
    - Linux: https://drive.google.com/file/d/0B_52viSz8FLNVmk2RUJOeU82M2M/view?usp=sharing
    - Mac: https://drive.google.com/file/d/0B_52viSz8FLNc2VvTERsUWVyVUU/view?usp=sharing

2. Open the terminal application and unarchive the file:

    ```shell
    unzip Installer_Linux.install.zip
    ```

3. Install FocalScan, including the MATLAB runtime (it is important that the correct runtime version is used, in this case v901):

    ```shell
    ./Installer_Linux.install
    ```
    or

    ```shell
    open ./Installer_Mac.app
    ```
    This should bring up a window for downloading and installing the runtime. Make a note of where FocalScan is installed and where the runtime is installed (such as /Applications/FocalScan/application and ∼/bin/MCR/v901, respectively, or any other directories chosen). The path to the runtime directory is from now on referred to as “MCR path”. 
    
    **Note:** Remote installation on a Linux server was observed to occassionally fail if X forwarding was not used (a bug in the MATLAB installer program): If no graphical interface appears, a text might still be printed suggesting that installation has "Finished". However, it might still **not** have succeded, unless the user has write access to the directory /usr. It also appears that the installer fails to display the graphical interface if a certain amount of time has passed after logging in to a remote computer using X forwarding. In this case, try again after reconnecting to the remote machine. 

    **Note:** If the installer fails to start on Linux, try changing permissions by typing ```chmod 755 Installer_Linux.install```.

4. Setup environment variables (in the following, substitute path to focalscan for the directory where the files were downloaded (for instance, ∼/focalscan):

    ```shell
    export PATH=$PATH:path_to_focalscan
    ```

    Example:

    ```shell
    export PATH=$PATH:/Applications/FocalScan/application
    ```

    **Important:** The “**application**” subdirectory needs to be included in the path.

5. Download and unzip example data to test the installation with, available at the following links:

    **Annotation files:**

    - https://drive.google.com/file/d/0B_52viSz8FLNTk4xb2VuY3lSQk0/view?usp=sharing

    **Test data:**

    - https://drive.google.com/file/d/0B_52viSz8FLNX2lZVF9SOTFMQzg/view?usp=sharing

6. Test the installation (gene-level analysis on breast cancer data from TCGA). In the following command, substitute “MCR root” for the installation directory of the MATLAB compiler runtime (for instance /Applications/MATLAB/MATLAB_Runtime/v901, if that is were it was installed):

    ```shell
    focalscan MCR_path expr_csv example_data/BRCA_expr.csv seg_file example_data/BRCA.seg annot_file annotation/gencode17_symbols.bed reportdir reports
    ```

    This might take up to 40 minutes or slightly longer, depending on processor speed and available memory.

    If the following occurs, double check that the correct input file paths were given:
    ```shell
    Invalid file identifier.  Use fopen to generate a valid file identifier.
    ```

7. Inspect the output:

    ```shell
    cat reports/peaks.txt
    ```

Any platform, usage from within the MATLAB environment (the function "FocalScan")
------------------------------------------------------

1. Open the terminal application and download the scripts (assuming that git is installed, otherwise just download the zip file):

    ```shell
    git clone https://github.com/jowkar/focalscan.git
    ```

2. Download and unzip example data to test the installation with, available at the following links:

    **Annotation files:**

    - https://drive.google.com/file/d/0B_52viSz8FLNTk4xb2VuY3lSQk0/view?usp=sharing

    **Test data:**

    - https://drive.google.com/file/d/0B_52viSz8FLNX2lZVF9SOTFMQzg/view?usp=sharing

3. Open MATLAB and add the scripts to the MATLAB path (will have to repeated each time MATLAB is opened, unless the startup settings are also changed):

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
    FocalScan('expr_csv','example_data/BRCA_expr.csv','seg_file','example_data/BRCA.seg','annot_file','annotation/gencode17_symbols.bed','reportdir','reports')
    ```

    This might take up to 40 minutes or slightly longer, depending on processor speed and available memory. The results will be saved to the directory “reports”.

# 2. Input files

RNA-seq data should be provided as read counts per gene/tile. Gene-level read counts can be obtained with HTSeq (or equivalent) and tile-level counts with coverageBed (for instance using the included script "quantify_tiles.sh").

Two expression data input options exist:
----------------------------------------

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
focalscan.sh <parameter1> <parameter1_value> ... <parameterN> <parameterN_value>
```
Example:

Using CSV input:
```shell
focalscan.sh expr_csv example_data/BRCA_expr.csv seg_file example_data/BRCA.seg annot_file annotation/gencode17_symbols.bed
```

Using separate read count files for each sample:
```shell
focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .gene_counts annot_file annotation/gencode17.bed seg_file example_data/BRCA.seg
```
Tile-level analysis:
```shell
focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .tile_counts annot_file annotation/gencode17.bed seg_file example_data/BRCA.seg
```
- From within MATLAB:
```MATLAB
FocalScan('parameter1','parameter1_value',...,'parameterN','parameterN_value')
```
Example:
```MATLAB
FocalScan('expr_csv','example_data/BRCA_expr.csv','seg_file','example_data/BRCA.seg','annot_file','annotation/gencode17_symbols.bed')
```

- With the compiled executable (no MATLAB installation)

```shell
focalscan <MCR_path> <parameter1> <parameter1_value> ... <parameterN> <parameterN_value>
```

where <MCR_path> is the path to the MATLAB runtime (v901, downloaded by the included installer).

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
- "pearson_corr": pearson correlation coefficient

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
standalone_peakdetection.sh <MCR_path> example_data/test_CSV/report.txt annotation/gencode17_symbols.bed 0.7 fs_hp ./new_peaks.txt
```

Annotating a tile peak report:
------------------------------

By default, the program will not write the IDs of genes overlapping the tiles in the output report unless a separate gene annotation file was added with the parameter "optional_gene_annot". To add this information afterwards, the script "annotate_peaks.sh" can be used:

```shell
annotate_peaks.sh peak_file_path annot_file_path out_file
```

Example:

```shell
annotate_peaks.sh peaks.txt annotation/gencode17_symbols.bed peaks_annotated.txt
```

The corresponding compiled version is named "annotate_peaks" (without .sh).
