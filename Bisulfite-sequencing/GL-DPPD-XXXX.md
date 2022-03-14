# Bioinformatics pipeline for bisulfite sequencing (methylseq) data

> **This document holds an overview and some example commands of how GeneLab processes bisulfite sequencing (methylseq) datasets. Exact processing commands for specific datasets that have been released are available in this repository [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory and are also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:**   
**Revision:** -  
**Document Number:** GL-DPPD-XXXX  

**Submitted by:**  
Michael D. Lee (GeneLab Analysis Team)  

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Jonathan Galazka (GeneLab Project Scientist)  

---

# Table of contents

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [1. Raw Data QC](#1-raw-data-qc)
  - [2. Adapter trimming/quality filtering](#2-adapter-trimmingquality-filtering)
    - [If not RRBS or if RRBS using MseI digestion](#if-not-rrbs-or-if-rrbs-using-msei-digestion)
    - [If RRBS with MspI digestion](#if-rrbs-with-mspi-digestion)
    - [If RRBS with NuGEN ovation kit](#if-rrbs-with-nugen-ovation-kit)
  - [3. Filtered/Trimmed Data QC](#3-filteredtrimmed-data-qc)
  - [4. Alignment](#4-alignment)
    - [Generate reference](#generate-reference)
    - [Align](#align)
  - [5. Deduplicate (only if not RRBS data)](#5-deduplicate-if-not-rrbs-data)
  - [6. Extract methylation calls](#6-extract-methylation-calls)
  - [7. Generate individual sample report](#7-generate-individual-sample-report)
  - [8. Generate combined summary report](#8-generate-combined-summary-report)
  - [9. Alignment QC](#alignment-qc)
  - [10. Generate MultiQC project report](#10-generate-multiqc-project-report)

---

# Software used

|Program|Version|Relevant Links|
|:------|:-----:|------:|
|FastQC| 0.11.9 |[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC| 1.10.1 |[https://multiqc.info/](https://multiqc.info/)|
|TrimGalore!| 0.6.6 |[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|Bismark| 0.23.0 |[https://github.com/FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark)|
|bowtie2| 2.4.2 |[https://github.com/BenLangmead/bowtie2#overview](https://github.com/BenLangmead/bowtie2#overview)|
|samtools| 1.11 |[https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)|
|qualimap| 2.2.2d |[http://qualimap.conesalab.org/](http://qualimap.conesalab.org/)|

---

# General processing overview with example commands

> Exact processing commands for specific datasets are available in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory of this repository, as well as being provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).  

## 1. Raw Data QC

```
fastqc -o raw_fastqc_output *raw.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*raw.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input data:**

* *raw.fastq.gz (raw reads)

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)


### 1a. Compile Raw Data QC

```
multiqc -o raw_multiqc_output -n raw_multiqc -z raw_fastqc_output/
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`-n` – the filename prefix of results
*	`-z` – specifies to zip the output data directory
*	`raw_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* raw_fastqc_output/*fastqc.zip (FastQC output data)

**Output data:**

* raw_multiqc_output/raw_multiqc_report.html (multiqc output html summary)
* raw_multiqc_output/raw_multiqc_data.zip (zipped directory containing multiqc output data)

<br>  

---

## 2. Adapter trimming/quality filtering
See `trim_galore --help` [menu](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3035) for more info on any of the below.

### If not RRBS or if RRBS using MseI digestion
Note that the `--rrbs` option is **not** appropriate when RRBS (reduced representation bisulfite sequencing) libraries were prepared with MseI digestion (see `trim_galore --help` menu [(starting at this line)](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3337).

**Single-end example**  

```bash
trim_galore --cores 4 --gzip sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip --cores 4 --paired sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw_trimmed.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_trimmed.fq.gz sample-1_R2_trimmed.fastq.gz
```

### If RRBS with MspI digestion
Note that if the library preparation was non-directional, we need to also add `--non_directional` to this command (whether single-end or paired-end). 

**Single-end example**  

```bash
trim_galore --gzip --cores 4 --rrbs sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip --cores 4 --rrbs --paired sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw_trimmed.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_trimmed.fq.gz sample-1_R2_trimmed.fastq.gz
```


### If RRBS with NuGEN ovation kit
Libraries prepared with the NuGEN ovation kit need to be procesed with an additional script provided by the company's [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq). 

Following their instructions, we first run an adapter-trimming/quality-filtering step with trimgalore. Note that the `--rrbs` option is not appropriate to pass to trimgalore when this kit is used (see `trim_galore --help` menu [(starting at this line)](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3329). Then we utilize the company's script to remove the random diversity sequences added by the kit. 

#### First adapter-trimming/quality-filtering with trimgalore

**Single-end example**  

```bash
trim_galore --cores 4 -a AGATCGGAAGAGC sample-1_raw.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --cores 4 --paired -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz
```

#### Now running NuGEN-specific script

We can download the script from their [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq) with the following:

```bash
curl -LO https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/trimRRBSdiversityAdaptCustomers.py
```

**Single-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py -1 sample-1_raw_trimmed.fq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_raw_trimmed.fq_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py -1 sample-1_R1_raw.fastq.gz -2 sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw.fastq_trimmed.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw.fastq_trimmed.fq.gz sample-1_R2_trimmed.fastq.gz
```

**Parameter Definitions for `trim_galore`:**  

* `--gzip` - gzip compress the output(s)
* `--cores` - number of cores to use
* `--rrbs` - specific trimming suitable for RRBS data generated with MspI digestion only
* `-a` - specific adapter sequence to be trimmed off of forward or single reads
* `-a2` - specific adapter sequence to be trimmed off of reverse reads
* `--paired` - specifies data are paired-end
* positional arguments represent the input files, 2 of them if paired-end data


**Parameter Definitions for `trimRRBSdiversityAdaptCustomers.py `:**  

- `-1` - forward or single input reads
- `-2` - reverse reads if paired-end data

**Input Data:**

* gzip-compressed fastq files (original reads)

**Output Data:**

* gzip-compressed fastq files (adapter-trimmed/quality-filtered reads)

<br>

---

## 3. Filtered/Trimmed Data QC
```
fastqc -o trimmed_fastqc_output/ *trimmed.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`*trimmed.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input data:**

* *trimmed.fastq.gz (filtered/trimmed reads)

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)

### 3a. Compile Trimmed Data QC

```
multiqc -o trimmed_multiqc_output -n trimmed_multiqc -z trimmed_fastqc_output/
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`-n` – the filename prefix of results
*	`-z` – specifies to zip the output data directory
*	`trimmed_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* trimmed_fastqc_output/*fastqc.zip (FastQC output data)

**Output data:**

* trimmed_multiqc_output/trimmed_multiqc_report.html (multiqc output html summary)
* trimmed_multiqc_output/trimmed_multiqc_data.zip (zipped directory containing multiqc output data)

<br>

---

## 4. Alignment

### Generate reference


### Align

<br>

---

## 5. Deduplicate (only if not RRBS data)

<br>

---

## 6. Extract methylation calls

<br>

---

## 7. Generate individual sample report

<br>

---

## 8. Generate combined summary report

<br>

---

## 9. Alignment QC

<br>

---

## 10. Generate MultiQC project report




---
