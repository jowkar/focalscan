# About

FocalScan identifies genomic regions where many tumors show simultaneous in- creases in DNA copy-number amplitude (CNA) and RNA expression (or conversely for DNA deletions). Empirically, many important oncogenes show this pattern of alteration. The FocalScan score is based on the dot product of CNA and RNA changes. This puts equal weight to the two variables, but requires coordinated changes in both to achieve a positive score.

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

- MATLAB

Recommeded for pre-processing expression data:

- bedtools (http://bedtools.readthedocs.org/en/latest/)
- samtools (http://www.htslib.org)
- HTseq (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)

Recommended for visualizing the results:

- IGV (https://www.broadinstitute.org/igv/)


# 1. Installation

1) Download the files. Example data can be found at:
- Annotation files: https://transfer.sh/8EZ7g/annotation.zip 
- Gene expression and copy number data: https://transfer.sh/nS0Qw/example-data.zip

To download the above files via the command line, use for instance (on Linux/Mac):
```curl https://transfer.sh/8EZ7g/annotation.zip -o annotation.zip```

2) Add execution permissions to the shell scripts by entering the directory and typing ```chmod +x *.sh```. (On Unix/Linux).

3) Either run the program from the same directory, or add the directory to the path (```export /path/to/FocalScan``` (Unix/Linux) or ```addpath(genpath('path/to/FocalScan'))``` (from within MATLAB)).

# 2. Input files

RNA-seq data should be provided as read counts per gene/tile. Gene-level read counts can be obtained with HTSeq (or equivalent) and tile-level counts with coverageBed (for instance using the included script "quantify_tiles.sh").

Two FocalScan expression data input options exist:
--------------------------------------------------

1) A CSV-file with samples as columns and genes as columns.

2) A combination of the following:
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
```
./quantify_tiles.sh <BAM_file> annotation/hg18_hg19_1kb_tiles.bed
```
This needs to be done for each .bam file, preferably in parallell to speed up processing.
A file with RNA-seq read counts (\*.tile_counts) for each genomic tile will be generated.

NOTE: Your .bam files may use chromosome names formatted as e.g. 'chr1' or simply '1'. In
the latter case, instead use the hg18_hg19_1kb_tiles_nochr.bed file.

NOTE: You may consider pre-filtering your .bam files to only consider uniquely mapped/high
quality reads (e.g. quality 255 only for TopHat alignments, by running 'samtools view -b
-q255 in.bam > out.bam').

NOTE: FocalScan does not take RNA-seq strand information into account. E.g. TCGA RNA-seq
datasets are not strand specific, but this could be useful in other cases. Strand-specific
analysis can be accomplished by first splitting .bam files into '+' and '-'
fractions using samtools, and running FocalScan on each fraction.

## 2) Run FocalScan:
-----------------

- Using the shell script:
```
./focalscan.sh <parameter1> <parameter1_value> ... <parameterN> <parameterN_value>
```
Example:

Using CSV input:
```
./focalscan.sh expr_csv example_data/BRCA_expr.csv seg_file example_data/BRCA_cna.seg annot_file annotation/gencode17_symbols.bed
```

Using separate read count files for each sample:
```
./focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .gene_counts annot_file annotation/gencode17.bed seg_file example_data/BRCA_cna.seg
```
Tile-level analysis:
```
./focalscan.sh expr_path example_data/read_count_files index_file example_data/index.txt file_extension .tile_counts annot_file annotation/gencode17.bed seg_file example_data/BRCA_cna.seg
```
- From within MATLAB:
```
FocalScan('parameter1','parameter1_value',...,'parameterN','parameterN_value')
```
Example:
```
FocalScan('expr_csv','example_data/BRCA_expr.csv','seg_file','example_data/BRCA_cna.seg','annot_file','annotation/gencode17_symbols.bed')
```
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

```./standalone_peakdetection.sh report_file_path annot_file_path peak_level scorefield out_file```

Valid options for the "scorefield" parameter (the metric to use as basis for peak detection) are:

- "fs_hp": the standard FocalScan score (with focality filter)
- "fs": FocalScan score without focality filter
- "sum_cna_hp": summed copy number amplitudes, with focality filter
- "sum_cna": summed copy number amplitudes, without focality filter
- "spearman_corr": spearman correlation coefficient

(Assuming that all of the above scores are present in the report.txt file, as is the case by default.)

Example:

```./standalone_peakdetection.sh example_data/test_CSV/report.txt annotation/gencode17_symbols.bed 0.7 fs_hp ./new_peaks.txt```

Or the MATLAB function:

```standalone_peakdetection(report_file_path,annot_file_path,peak_level,scorefield,out_file)```
