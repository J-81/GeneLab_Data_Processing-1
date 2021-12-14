# GeneLab raw data generation for bulk RNAseq prepared with the Qiagen UPX kit

> **This page holds an overview and instructions for how GeneLab generates raw RNA sequence data from libraries prepared using the [Qiagen UPX kit](https://www.qiagen.com/us/products/discovery-and-translational-research/next-generation-sequencing/rna-sequencing/three-rnaseq/qiaseq-upx-3-transcriptome-kits/). The instructions below assume that the samples were prepared as follows:
> 1. Bulk RNA was was extracted from samples 
> 2. Each sample was prepared with a unique cell ID
> 3. Samples were pooled together to create sample pools
> 4. Each sample pool was prepared with a unique single index
> 5. Sample pools were pooled together to create a pool of pooled samples
> 6. The pool of pooled samples was sequenced on an Illumina NovaSeq instrument**  

---

**Date:** September 13, 2021  
**Revision:** D  
**Document Number:** GL-DPPD-7101-D  

**Submitted by:**  
Amanda Saravia-Butler (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)   
Jonathan Galazka (GeneLab Project Scientist)

---

## Updates from previous revision

Two additional sub-steps were added to step 4:
- A step to compile the alignment log files using multiQC, [step 4b](#4b-compile-alignment-logs)
- A step to index the alignment files, [step 4c](#4c-index-aligned-reads), which is required to assess read strandedness

Two additional steps were added prior to aligned read quantitation:
- Step 5, [5a](#5a-convert-gtf-to-genepred-file) and [5b](#5b-convert-genepred-to-bed-file), was added to create a reference annotation BED file required to assess read strandedness
- Step 6 was added to [determine read strandedness](#6a-determine-read-strandedness) and [compile read strandedness reports](#6b-compile-strandedness-reports), to determine which RSEM `--strandedness` setting to use during aligned read quantitation

The aligned read quantitation step, now [step 8](#8-count-aligned-reads-with-rsem), was modified to use the results of the read strandedness step to inform the correct RSEM `--strandedness` setting
> Note: A subset of samples from all datasets previously processed were evaluated for strandedness, and those datasets identified to have been processed with the incorrect RSEM `--strandedness` setting were reprocessed with the correct setting

The DESeq2 Normalization and DGE step for datasets with ERCC spike-in, now [step 9a](#9a-for-datasets-with-ercc-spike-in), was modified as follows:
- Perform ERCC normalization using only ERCC group B genes, since the concentration of these genes are the same in ERCC mix 1 and mix 2
- Remove any samples that do not contain detectible ERCC group B spike-ins prior to generation and subsequent analysis of ERCC-normalized count data
- Account for the edge case in which rescaling using ERCC size factors fails due to zero gene counts

Added "Stat_" column containing the Wald Statistic (similar to a Z-score) to the DGE output tables for datasets both with and without ERCC spike-in, now [step 9a](#9a-for-datasets-with-ercc-spike-in) and [step 9b](#9b-for-datasets-without-ercc-spike-in), respectively, which will be used for GSEA visualizations

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - **1. Raw Data QC**
    - [**1a. Raw Data QC**](#1a-raw-data-qc)
    - [**1b. Compile Raw Data QC**](#1b-compile-raw-data-qc)
  - **2. Trim/Filter Raw Data and Trimmed Data QC**
    - [**2a. Trim/Filter Raw Data**](#2a-trimfilter-raw-data)
    - [**2b. Trimmed Data QC**](#2b-trimmed-data-qc)
    - [**2c. Compile Trimmed Data QC**](#2c-compile-trimmed-data-qc)
  - [**3. Build STAR Reference**](#3-build-star-reference)
  - **4. Align Reads to Reference Genome then Index**
    - [**4a. Align Reads to Reference Genome with STAR**](#4a-align-reads-to-reference-genome-with-star)
    - [**4b. Compile Alignment Logs**](#4b-compile-alignment-logs)
    - [**4c. Index Aligned Reads**](#4c-index-aligned-reads)
  - **5. Create Reference BED File**
    - [**5a. Convert GTF to genePred File**](#5a-convert-gtf-to-genepred-file)
    - [**5b. Convert genePred to BED File**](#5b-convert-genepred-to-bed-file)
  - **6. Assess Strandedness**
    - [**6a. Determine Read Strandedness**](#6a-determine-read-strandedness)
    - [**6b. Compile Strandedness Reports**](#6b-compile-strandedness-reports)
  - [**7. Build RSEM Reference**](#7-build-rsem-reference)
  - [**8. Count Aligned Reads with RSEM**](#8-count-aligned-reads-with-rsem)
  - [**9. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R**](#9-normalize-read-counts-perform-differential-gene-expression-analysis-and-add-gene-annotations-in-r)
    - [**9a. For Datasets with ERCC Spike-In**](#9a-for-datasets-with-ercc-spike-in)
    - [**9b. For Datasets without ERCC Spike-In**](#9b-for-datasets-without-ercc-spike-in)
  
---

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|FastQC|`fastqc -v`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`multiqc -v`|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|`cutadapt --version`|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|`trim_galore -v`|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|`STAR --version`|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|RSEM|`rsem-calculate-expression --version`|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Samtools|`samtools --version`|[http://www.htslib.org/](http://www.htslib.org/)|
|gtfToGenePred|N/A|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed|N/A|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|infer_experiment|`infer_experiment.py --version`|[http://rseqc.sourceforge.net/#infer-experiment-py](http://rseqc.sourceforge.net/#infer-experiment-py)|
|Bioconductor|`BiocManager::version()`|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|`packageVersion("DESeq2")`|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|`packageVersion("tximport")`|[https://bioconductor.org/packages/release/bioc/html/tximport.html](https://bioconductor.org/packages/release/bioc/html/tximport.html)|
|tidyverse|`packageVersion("tidyverse")`|[https://www.tidyverse.org](https://www.tidyverse.org)|
|Risa|`packageVersion("Risa")`|[https://www.bioconductor.org/packages/release/bioc/html/Risa.html](https://www.bioconductor.org/packages/release/bioc/html/Risa.html)|
|STRINGdb|`packageVersion("STRINGdb")`|[https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|PANTHER.db|`packageVersion("PANTHER.db")`|[https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html)|
|org.Hs.eg.db|`packageVersion("org.Hs.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)|
|org.Mm.eg.db|`packageVersion("org.Mm.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html)|
|org.Dm.eg.db|`packageVersion("org.Dm.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html)|
|org.Ce.eg.db|`packageVersion("org.Ce.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html)|
|org.At.tair.db|`packageVersion("org.At.tair.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)|
|org.EcK12.eg.db|`packageVersion("org.EcK12.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html)|
|org.Sc.sgd.db|`packageVersion("org.Sc.sgd.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)|

>**\*** Exact versions are available along with the processing commands for each specific dataset in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory. 

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are provided in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory and are also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).  

---

## 1a. Raw Data QC  

```
fastqc -o /path/to/raw_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**
- *fastq.gz (raw reads)

**Output Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

## 1b. Compile Raw Data QC  

```
multiqc -n raw_multiqc -o /path/to/raw_multiqc/output/directory /path/to/directory/containing/raw_fastqc/files
```

**Parameter Definitions:**

* `-n` - prefix name for output files
* `-o` – the output directory to store results
* `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

**Output Data:**
- raw_multiqc.html (multiqc report)
- raw_multiqc_data (directory containing multiqc data)

<br>

---

## 2a. Trim/Filter Raw Data  

```
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --phred33 \
  --illumina \ # if adapters are not illumina, replace with adapters used
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \ # only for PE studies, remove this paramater if raw data are SE
  sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz
# if SE, replace the last line with only the forward reads (R1) of each sample

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
