# GeneLab raw data generation for bulk RNAseq prepared with the Qiagen UPX kit

> **This page holds an overview and instructions for how GeneLab generates raw RNA sequence data from libraries prepared using the [Qiagen UPX kit](https://www.qiagen.com/us/products/discovery-and-translational-research/next-generation-sequencing/rna-sequencing/three-rnaseq/qiaseq-upx-3-transcriptome-kits/). The instructions below assume that the samples were prepared as follows:**
> 1. Bulk RNA was was extracted from samples 
> 2. Each sample was prepared with a unique cell ID
> 3. Samples were pooled together to create sample pools
> 4. Each sample pool was prepared with a unique single index
> 5. Sample pools were pooled together to create a pool of pooled samples
> 6. The pool of pooled samples was sequenced on an Illumina NovaSeq instrument

---

**Date:** September 3, 2021  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX-  

**Submitted by:**  
Amanda Saravia-Butler (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)   
Jonathan Galazka (GeneLab Project Scientist)

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - **1. Demultiplex Sample Pools**
    - [**1a. Raw Data QC**](#1a-raw-data-qc)
    - [**1b. Compile Raw Data QC**](#1b-compile-raw-data-qc)
  - **2. Identify Cell IDs**
    - [**2a. Trim/Filter Raw Data**](#2a-trimfilter-raw-data)
    - [**2b. Trimmed Data QC**](#2b-trimmed-data-qc)
    - [**2c. Compile Trimmed Data QC**](#2c-compile-trimmed-data-qc)
  - [**3. Extract Cell IDs and UMIs**](#3-build-star-reference)
  - **4. Demultiplex Individual Samples**
    - [**4a. Align Reads to Reference Genome with STAR**](#4a-align-reads-to-reference-genome-with-star)
    - [**4b. Compile Alignment Logs**](#4b-compile-alignment-logs)
    - [**4c. Index Aligned Reads**](#4c-index-aligned-reads)
  
---

# Software used  

|Program|Version|Relevant Links|
|:------|:-----:|:-------------|
|bcl2fastq|2.20|[https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)|
|umi_tools|1.1.2|[https://umi-tools.readthedocs.io/en/latest/](https://umi-tools.readthedocs.io/en/latest/)|

---

# General processing overview with example commands   

---

## 1. Demultiplex Sample Pools  

```
bcl2fastq --runfolder-dir /path/to/NovaSeq/directory \
  --output-dir /path/to/output/directory \
  --loading-threads 8 \
  --writing-threads 8 \
  --processing-threads 12 \
  --barcode-mismatches 1 \
  --no-lane-splitting \
  --no-bgzf-compression \
  --minimum-trimmed-read-length 0 \
  --mask-short-adapter-reads 0
```

**Parameter Definitions:**

* `--runfolder-dir` – path to the directory output from the NovaSeq sequencing run  
* `--output-dir` – path to the directory to store the bcl2fastq output files 
* `--loading-threads` - number of threads used to load the bcl data
* `--writing-threads` - number of threads used to write the fastq data
* `--processing-threads` - number of threads used to demultiplex the data
* `--barcode-mismatches` - number of mismatches allowed per index adapter
* `--no-lane-splitting` - instructs the program not to split the fastq files by lane (only use if the same pool of pooled samples is run on all lanes of the flow cell)
* `--no-bgzf-compression` - turn off bgzf and instead use gzip to compress the fastq files
* `--minimum-trimmed-read-length` - specifies the minimum read length after adapter trimming, a value of 0 indicates no adapter trimming
* `--mask-short-adapter-reads` - specifies the number of reads to mask if less than the minimum trimmed read length

**Input Data:**
- NovaSeq/directory (output directory from the NovaSeq run - must contain the SampleSheet.csv file)
   > Note: Make sure the NovaSeq output directory contains a properly formatted [SampleSheet.csv](https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/sequencing-sheet-format-specifications-technical-note-970-2017-004.pdf) file in which the samples listed are the sample pools rather than individual samples.  

**Output Data:**
- *fastq.gz (fastq.gz files for each sample pool listed in the SampleSheet.csv file plus a set of fastq.gz files for Undetermined reads)
- /Stats (directory containing demultiplexing statistics)
- /Reports (directory containing html reports of the demultiplexing statistics)

<br>

## 2. Identify Cell IDs  

```
umi_tools whitelist --stdin=/path/to/${sample_pool}*R2*fastq.gz \
	--bc-pattern=CCCCCCCCCCNNNNNNNNNNNN \
	--plot-prefix=/path/to/plot/output/files/${sample_pool} \
	-L /path/to/${sample_pool}_whitelist.log \
	--set-cell-number Number_of_cellIDs \
	--subset-reads Number_of_reads > /path/to/whitelist/output/files/${sample_pool}_whitelist.tsv
```

**Parameter Definitions:**

* `--stdin` - read contining the cell ID (this is the reverse read for samples prepared with the Qiagen UPX kit) 
* `--bc-pattern` – pattern for the barcode on the read containing the cell ID; `C`s indicate placeholders for cell IDs; `N`s indicate placeholders for UMIs 
* `--plot-prefix` - instructs the program to output plots to visualise the set of thresholds considered for defining cell barcodes
* `-L` - specifies whitelist output log file
* `--set-cell-number` - number of cell IDs to extract; this should match the number of individual samples in each sample pool
* `--subset-reads` - number of reads to use to identify true cell barcodes; to use all reads set this number to greater than the max number of reads
* `${sample_pool}_whitelist.tsv` - specifies the file to output the cell IDs identified for each sample pool

**Input Data:**
- *R2\*fastq.gz (reverse fastq.gz file for each sample pool, generated from step 1)

**Output Data:**
- *whitelist.log (whitelist extraction log file)
- *barcode_counts.png (plot visualizing the knee in the cell barcode count distrubution)
- *barcode_knee.png (plot visualizing the knee in the cell barcode distrubution)
- *whitelist.tsv (whitelist output file containing 4 tab-separated columns: 1-whitelist cellID, 2-cellIDs that are 1bp different from the respective whitelist cellID identified, 3-number of whitelisted cellIDs, 4-number of 1bp different cellIDs)
  > **Note1:** Review column 1 of this file for each sample pool and check that all cellIDs for each sample in the respective sample pool are present. If any cellIDs are missing, increase the `--set-cell-number` option by one, and re-run this step. Repeat until all cellIDs for each sample pool are shown in column 1. 
  >
  > **Note2:** Before moving on to the next step, if any cellIDs are present in column 1 that are not associated with a sample in the sample pool, remove that row(s) from the *whitelist.tsv file before moving on to the next step (it's good proactice to keep a copy of the original *whitelist.tsv before modifying the file).

<br>

---

## 3. Extract Cell IDs and UMIs  

```
umi_tools extract --stdin=./Fastq/${sample}_R2_001.fastq.gz \
	--read2-in=./Fastq/${sample}_R1_001.fastq.gz \
	--bc-pattern=CCCCCCCCCCNNNNNNNNNNNN \
	--stdout=./Fastq_cellID-umi-extracted/${sample}_R1_raw.fastq.gz \
	--read2-stdout \
	--whitelist=./whitelist_files/${sample}_whitelist.tsv \
	--error-correct-cell \
	--filter-cell-barcode \
	--log=./cellID-umi-extraction_logs/${sample}_R1_raw_extraction.log
```

**Parameter Definitions:**

* `--gzip` – compress the output files with `gzip`
* `--path_to_cutadapt` - specify path to cutadapt software if it is not in your `$PATH`
* `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
* `--illumina` - defines the adapter sequence to be trimmed as the first 13bp of the Illumina universal adapter `AGATCGGAAGAGC`
* `--output_dir` - the output directory to store results
* `--paired` - indicates paired-end reads - both reads, forward (R1) and reverse (R2) must pass length threshold or else both reads are removed
* `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

**Input Data:**
- *fastq.gz (raw reads)

**Output Data:**
- *fastq.gz (trimmed reads)
- *trimming_report.txt (trimming report)

<br>

## 2b. Trimmed Data QC  

```
fastqc -o /path/to/trimmed_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**
- *fastq.gz (trimmed reads)

**Output Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

## 2c. Compile Trimmed Data QC  

```
multiqc -n trimmed_multiqc -o /path/to/trimmed_multiqc/output/directory /path/to/directory/containing/trimmed_fastqc/files
```

**Parameter Definitions:**

* `-n` - prefix name for output files
* `-o` – the output directory to store results
* `/path/to/directory/containing/trimmed_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

**Output Data:**
- trimmed_multiqc.html (multiqc report)
- trimmed_multiqc_data (directory containing multiqc data)

<br>

---

## 3. Build STAR Reference  
