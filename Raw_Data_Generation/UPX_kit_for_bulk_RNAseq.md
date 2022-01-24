# GeneLab Raw Data Generation For Bulk RNAseq Prepared With The Qiagen UPX Kit

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
  - **1. [Demultiplex Sample Pools](#1-demultiplex-sample-pools)**
  - **2. [Identify Cell IDs](#2-identify-cell-ids)**
  - **3. [Extract Cell IDs and UMIs](#3-extract-cell-ids-and-umis)**
  - **4. [Demultiplex Individual Samples](#4-demultiplex-individual-samples)**
  
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
	--loading-threads NumberOfThreads \
	--writing-threads NumberOfThreads \
	--processing-threads NumberOfThreads \
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
	--set-cell-number NumberOfCellIDs \
	--subset-reads NumberOfReads > /path/to/whitelist/output/files/${sample_pool}_whitelist.tsv
```

**Parameter Definitions:**

* `--stdin` - fastq file with the reads contining the cell ID and UMI (this is the reverse reads for samples prepared with the Qiagen UPX kit) 
* `--bc-pattern` – pattern for the barcode on the read containing the cell ID and UMI; `C`s indicate placeholders for cell IDs; `N`s indicate placeholders for UMIs 
* `--plot-prefix` - instructs the program to output plots to visualise the set of thresholds considered for defining cell barcodes
* `-L` - specifies whitelist output log file
* `--set-cell-number` - number of cell IDs to extract; this should match the number of individual samples in each sample pool
* `--subset-reads` - number of reads to use to identify true cell barcodes; to use all reads set this number to greater than the max number of reads
* `${sample_pool}_whitelist.tsv` - specifies the file to output the cell IDs identified for each sample pool

**Input Data:**
- *R2\*fastq.gz (reverse fastq.gz file for each sample pool, generated from [Step 1](#1-demultiplex-sample-pools))

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
umi_tools extract --stdin=/path/to/${sample_pool}*R2*fastq.gz \
	--read2-in=/path/to/${sample_pool}*R1*fastq.gz \
	--bc-pattern=CCCCCCCCCCNNNNNNNNNNNN \
	--stdout=/path/to/Fastq/output/files/${sample_pool}_R1_raw.fastq.gz \
	--read2-stdout \
	--whitelist=/path/to/whitelist/files/${sample_pool}_whitelist.tsv \
	--error-correct-cell \
	--filter-cell-barcode \
	--log=/path/to/extract/log/files/${sample_pool}_R1_raw_extraction.log
```

**Parameter Definitions:**

* `--stdin` - fastq file with the reads contining the cell ID and UMI (this is the reverse read for samples prepared with the Qiagen UPX kit) 
* `--read2-in` - fastq file with the reads contining the sequence of interest (this is the forward read for samples prepared with the Qiagen UPX kit)
* `--bc-pattern` – pattern for the barcode on the read containing the cell ID and UMI; `C`s indicate placeholders for cell IDs; `N`s indicate placeholders for UMIs 
* `--stdout` – specifies the path and file name of the output fastq file 
* `--read2-stdout` - instructs the program to only output the fastq file designated with the `--read2-in` option 
* `--whitelist` - specifies the *whitelist.tsv file for each sample pool 
* `--error-correct-cell` - instructs the program to correct any single basepair mismatches in the cellID identified in column 2 of the *whitlist.tsv file for each sample pool
* `--filter-cell-barcode` - instructs the program to filter cell barcodes according to those provided in the *whitelist.tsv file
* `--log` - specifies the file to output the umi_tools extraction logs


**Input Data:**
- *R2\*fastq.gz (reverse fastq.gz file for each sample pool, generated from [Step 1](#1-demultiplex-sample-pools))
- *R1\*fastq.gz (forward fastq.gz file for each sample pool, generated from [Step 1](#1-demultiplex-sample-pools))
- *whitelist.tsv (whitelist of accepted cell barcodes for each sample pool, output from [Step 2](#2-identify-cell-ids))

**Output Data:**
- *R1_raw.fastq.gz (output fastq file containing the reads of interest with the cell ID and UMI in the read header)
- *R1_raw_extraction.log (umi_tools extraction log file)

<br>

## 4. Demultiplex Individual Samples 

Create a variable for each cellID in each sample pool:
```
first=$(cat ./cellIDs/${sample_pool}_cellIDs.txt | sed -n 1p)
second=$(cat ./cellIDs/${sample_pool}_cellIDs.txt | sed -n 2p)
third=$(cat ./cellIDs/${sample_pool}_cellIDs.txt | sed -n 3p)
...
```

Check that each cellID is defined correctly by the respective variable:
```
echo "First_cellID: ${first}"
echo "Second_cellID: ${second}"
echo "Third_cellID: ${third}"
...
```

Use the cellID in each read header to parse the reads in the sample pool fastq file to create an individual fastq file for each sample:
```
zcat /path/to/sample_pool/fastq/files/${sample_pool}_R1_raw.fastq.gz | sed -n '/_${first}_/{p; n;p; n;p; n;p}' | gzip > /path/to/individual/sample/fastq/output/files/${sample_pool}_1_${first}_R1_raw.fastq.gz

zcat /path/to/sample_pool/fastq/files/${sample_pool}_R1_raw.fastq.gz | sed -n '/_${second}_/{p; n;p; n;p; n;p}' | gzip > /path/to/individual/sample/fastq/output/files/${sample_pool}_2_${second}_R1_raw.fastq.gz

zcat /path/to/sample_pool/fastq/files/${sample_pool}_R1_raw.fastq.gz | sed -n '/_${third}_/{p; n;p; n;p; n;p}' | gzip > /path/to/individual/sample/fastq/output/files/${sample_pool}_3_${third}_R1_raw.fastq.gz

...

```

**Input Data:**
- *cellIDs.txt (single column list of each cellID in the respective sample pool)
- *R1_raw.fastq.gz (sample pool fastq file containing the reads of interest with the cell ID and UMI in the read header, output from [Step 3](#3-extract-cell-ids-and-umis))
- ... = repeat for all cellIDs in each sample pool

**Output Data:**
- *cellID_R1_raw.fastq.gz (fastq file containing reads from an indicvidual sample within a sample pool)
  > **Note:** After all sample pool fastq files have been parsed, individual sample fastq files can be renamed  

<br>
