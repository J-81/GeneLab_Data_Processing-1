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

- [Bioinformatics pipeline for bisulfite sequencing (methylseq) data](#bioinformatics-pipeline-for-bisulfite-sequencing-methylseq-data)
- [Table of contents](#table-of-contents)
- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Raw Data QC](#1-raw-data-qc)
    - [1a. Compile Raw Data QC](#1a-compile-raw-data-qc)
  - [2. Adapter trimming/quality filtering](#2-adapter-trimmingquality-filtering)
    - [If not RRBS or if RRBS using MseI digestion](#if-not-rrbs-or-if-rrbs-using-msei-digestion)
    - [If RRBS with MspI digestion](#if-rrbs-with-mspi-digestion)
    - [If RRBS with NuGEN ovation kit](#if-rrbs-with-nugen-ovation-kit)
      - [First adapter-trimming/quality-filtering with trimgalore](#first-adapter-trimmingquality-filtering-with-trimgalore)
      - [Now running NuGEN-specific script](#now-running-nugen-specific-script)
  - [3. Filtered/Trimmed Data QC](#3-filteredtrimmed-data-qc)
    - [3a. Compile Trimmed Data QC](#3a-compile-trimmed-data-qc)
  - [4. Alignment](#4-alignment)
    - [4a. Generate reference](#4a-generate-reference)
    - [4b. Align](#4b-align)
  - [5. Deduplicate (skip if data are RRBS)](#5-deduplicate-skip-if-data-are-rrbs)
  - [6. Extract methylation calls](#6-extract-methylation-calls)
  - [7. Generate individual sample report](#7-generate-individual-sample-report)
  - [8. Generate combined summary report](#8-generate-combined-summary-report)
  - [9. Alignment QC](#9-alignment-qc)
  - [10. Generate MultiQC project report](#10-generate-multiqc-project-report)
  - [11. Generate reference genome annotation information](#11-generate-reference-genome-annotation-information)
    - [11a. GTF to BED conversion](#11a-gtf-to-bed-conversion)
    - [11b. Making a mapping file of genes to transcripts](#11b-making-a-mapping-file-of-genes-to-transcripts)
    - [11c. Making a table of all gene annotations](#11c-making-a-table-of-all-gene-annotations)
  - [12. Differential methylation analysis](#12-differential-methylation-analysis)

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
|methylKit|1.20.0|[https://bioconductor.org/packages/release/bioc/html/methylKit.html](https://bioconductor.org/packages/release/bioc/html/methylKit.html)|

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
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

### If RRBS with MspI digestion
Note that if the library preparation was non-directional, the `--non_directional` flag needs to be added to this command (whether single-end or paired-end). 

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
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```


### If RRBS with NuGEN ovation kit
Libraries prepared with the NuGEN ovation kit need to be procesed with an additional script provided by the company's [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq). 

Following their instructions, we first run an adapter-trimming/quality-filtering step with trimgalore. Note that the `--rrbs` option is not appropriate to pass to trimgalore when this kit is used (see `trim_galore --help` menu [(starting at this line)](https://github.com/FelixKrueger/TrimGalore/blob/072ecf9a1f80f9eb41c8116c32284492f481cbbb/trim_galore#L3329). Then we utilize the company's script to remove the random diversity sequences added by the kit. 

#### First adapter-trimming/quality-filtering with trimgalore

**Single-end example**  

```bash
trim_galore --cores 4 -a AGATCGGAAGAGC sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --cores 4 --paired \
            -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC \
            sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

#### Now running NuGEN-specific script

The NuGEN-specific script can be downloaded from their [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq) with the following:

```bash
curl -LO https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/trimRRBSdiversityAdaptCustomers.py
```

**Single-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py -1 sample-1_trimmed.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_trimmed.fastq_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py \
        -1 sample-1_R1_trimmed.fastq.gz \
        -2 sample-1_R2_trimmed.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_trimmed.fastq_trimmed.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_trimmed.fastq_trimmed.fq.gz sample-1_R2_trimmed.fastq.gz
```

**Parameter Definitions for `trim_galore`:**  

* `--gzip` - gzip compress the output(s)
* `--cores` - number of cores to use
* `--rrbs` - specific trimming suitable for RRBS data generated with MspI digestion only
* `-a` - specific adapter sequence to be trimmed off of forward or single reads
* `-a2` - specific adapter sequence to be trimmed off of reverse reads
* `--paired` - specifies data are paired-end
* positional arguments represent the input read files, 2 of them if paired-end data


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

### 4a. Generate reference
The reference will need to be specific to what was sequenced. Bismark operates on a directory holding the target reference genome in fasta format.

```bash
# creating directory to hold reference and moving it into there
mkdir reference-genome
mv ref-genome.fasta reference-genome/

bismark_genome_preparation --parallel 4 reference-genome/
```

**Parameter Definitions:**

*	`--parallel` – specifies how many threads to use (note these will be doubled as it operates on both strands simultaneously)
*  positional argument specifing the directory holding the reference genome (should end in ".fa" or ".fasta", can be gzipped and including ".gz")

**Input data:**

* a directory holding the reference genome in fasta format

**Output data:**

* an additional subdirectory added to the the reference genome directory that was provided as input which holds indexes for the bisulfite converted reference genome

> **NOTE**  
> If using RNA, need to add the `--hisat` flag.

### 4b. Align

Note that if the library preparation was non-directional, the `--non_directional` flag needs to be added to this command (whether single-end or paired-end). 

**Single-end example**  

```bash
bismark --bam -p 4 --genome reference-genome/ sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
bismark --bam -p 4 --genome reference-genome/ \
        -1 sample-1_R1_trimmed.fastq.gz \
        -2 sample-1_R2_trimmed.fastq.gz
```

**Parameter Definitions:**

* `--bam` - specifies to convert the default output sam format into compressed bam format
* `-p` - allows us to specify the number of threads to use (will be doubled for operating on both strands simultaneously)
* `--genome` - specifies the directory holding the reference genome indexes (the same that was provided to the Generate reference step above)
* input trimmed-reads are provided as a positional argument if they are single-end data
* `-1` - where to specify the forward trimmed reads if paired-end data
* `-2` - where to specify the reverse trimmed reads if paired-end data


**Input data:**
* directory holding indexes of reference genome
* gzip-compressed fastq files (adapter-trimmed/quality-filtered reads)

**Output data:**  

* \*.bam - mapping file 
* \*_report.txt - bismark mapping report file


> **NOTE**  
> If using RNA, need to add the `--hisat` flag.


<br>

---

## 5. Deduplicate (skip if data are RRBS)
> **NOTE**  
> This step should **not** be done if the data are RRBS (reduced representation bisulfite sequencing; see [bismark documentation](https://github.com/FelixKrueger/Bismark/tree/master/Docs)).

```bash
deduplicate_bismark sample-1_trimmed_bismark_bt2.bam
```

**Parameter Definitions:**

* positional argument is the bam file produced in step 4 above

**Input data:**

* sample-1_trimmed_bismark_bt2.bam - mapping file produced in step 4 above

**Output data:**

* \*.deduplicated.bam - a deduplicated mapping file
* \*.deduplication_report.txt - report file of deduplication 


<br>

---

## 6. Extract methylation calls


**Single-end example**  

```bash
bismark_methylation_extractor --bedGraph --gzip --comprehensive sample-1_trimmed_bismark_bt2.bam
    # note, input should be the deduplicated version produced in step 5 above if not working with RRBS data
```

**Paired-end example**  

```bash
bismark_methylation_extractor --bedGraph --gzip --comprehensive --ignore_r2 2 --ignore_3prime_r2 2 sample-1_trimmed_bismark_bt2.bam
    # note, input should be the deduplicated version produced in step 5 above if not working with RRBS data
```


**Parameter Definitions:**

* `--bedGraph` - specifies to generate a bedGraph-formatted file of methylated CpGs (see bismark docs [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#optional-bedgraph-output)
* `--gzip` - specifies to gzip-compress teh larger output files
* `--comprehensive` - specifies to merge all four possible strand-specific methylation outputs into context-dependent output files
* `--ignore_r2` - allows specifying how many bases to ignore from the 5' end of the reverse reads (bismark docs recommend 2, see [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#appendix-iii-bismark-methylation-extractor))
* `--ignore_3prime_r2` - allows specifying how many bases to ignore from the 3' end of the reverse reads (this is utilized in the [nf-core methylseq workflow](https://nf-co.re/methylseq), set at [this line](https://github.com/nf-core/methylseq/blob/03972a686bedeb2920803cd575f4d671e9135af0/main.nf#L643)) 
* the positional argument is an input bam file

**Input data:**

* sample-1_trimmed_bismark_bt2.bam - bam file produced above (in step 4 if data are RRBS, or step 5 if not)

> **NOTE**  
> If data are **not** RRBS, the input bam file should be the deduplicated one produced by step 5 above. 


**Output data:**

* *\.txt.gz - bismark methylation call files for CpG, CHG, and CHH contexts (where H is A, T, or C) that were detected; see [bismark documentation](https://github.com/FelixKrueger/Bismark/tree/master/Docs), namely [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#methylation-call) for symbols, and [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#iv-bismark-methylation-extractor) for file format
* \*.bedGraph.gz - gzip-compressed bedGraph-formatted file of methylation percentages of each CpG site (see bismark docs [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#optional-bedgraph-output)
* \*.bismark.cov.gz - gzip-compressed bedGraph-formatted file like above "\*.bedGraph.gz", but also including 2 more columns of methylated and unmethylated counts at the specified position (see bismark docs [here](https://github.com/FelixKrueger/Bismark/tree/master/Docs#optional-bedgraph-output)
* \*_splitting_report.txt - text file of general detected methylation information
* \*.M-bias.txt - text file with methylation information in the context of the position in reads, helpful for investigating bias as a function of base position in the read

<br>

---

## 7. Generate individual sample report


```bash
bismark2report --alignment_report sample-1_trimmed_bismark_bt2_SE_report.txt \
               --splitting_report sample-1_trimmed_bismark_bt2_splitting_report.txt \
               --mbias_report sample-1_trimmed_bismark_bt2.M-bias.txt
```

**Parameter Definitions:**

* `--alignment_report` - where to provide the alignment report generated by step 4 above (e.g. "\*_bt2_SE_report.txt" if single-end data or "\*_bt2_PE_report.txt" if paired-end data)
* `--splitting_report` - where to provide the splitting report generated by step 6 above
* `--mbias_report` - where to provide the bias report generated by step 6 above

**Input data:**

* sample-1_trimmed_bismark_bt2_SE_report.txt - alignment report generated by step 4 above
* sample-1_trimmed_bismark_bt2_splitting_report.txt - splitting report generated by step 6 above
* sample-1_trimmed_bismark_bt2.M-bias.txt - bias report generated by step 6 above

> **NOTE**  
> If data are **not** RRBS, the deduplication report from step 5 above should also be provided to the above command, e.g.: `--dedup_report sample-1_trimmed_bismark_bt2.deduplication_report.txt` 

**Output data:**

* \*_report.html - a summary html file for the given sample

<br>

---

## 8. Generate combined summary report

```bash
bismark2summary 
```

**Input data:**  

* autodetects appropriate files in current working directory intially based on bam files generated in step 4 above

**Output data:**  

* bismark_summary_report.txt - summary table of general information on all samples
* bismark_summary_report.html - html summary of general information on all samples


<br>

---

## 9. Alignment QC

```bash
# sorting bam file
samtools sort -@ 4 -o sample-1_trimmed_bismark_bt2.sorted.bam \
         sample-1_trimmed_bismark_bt2.bam
    # note, input should be the deduplicated version produced 
    # in step 5 above if not working with RRBS data

qualimap bamqc -bam sample-1_trimmed_bismark_bt2.sorted.bam \
         -outdir sample-1_trimmed_bismark_bt2_qualimap \
         --collect-overlap-pairs --java-mem-size=6G -nt 4
```

**Parameter Definitions for `samtools`:**

* `sort` - specifies the sub-program of `samtools`
* `-@` - where to specify the number of threads to use
* `-o` - specifies the output file name
* the positional argument is the input bam file 

**Parameter Definitions for `qualimap`:**

* `bamqc` - specifies the sub-program of `qualimap`
* `-bam` - where to specify the input bam file
* `-outdir` - where to specify the output directory
* `--collect-overlap-pairs` - includes statistics of overlapping paired-end reads (if data were paired-end, no effect if single-end)
* `--java-mem-size=6G` - where to specify the amount of memory to use (here 6G)
* `-nt` - where to specify the number of threads to use

**Input data:**

* sample-1_trimmed_bismark_bt2.bam - bam file produced above (in step 4 if data are RRBS, or step 5 if not)

> **NOTE**  
> If data are **not** RRBS, the input bam file should be the deduplicated one produced by step 5 above. 

**Output data:**

* `sample-1_trimmed_bismark_bt2_qualimap/` - subdirectory of multiple alignment QC output files presented in an html file (see [qualimap documentation](http://qualimap.conesalab.org/doc_html/analysis.html#output))


<br>

---

## 10. Generate MultiQC project report

```bash
multiqc -o project_multiqc_output -n project_multiqc -z ./
```

**Parameter Definitions:**

*	`-o` – where to specify the output directory to store results
*	`-n` – where to specify the filename prefix of results
*	`-z` – specifies to zip the output data directory
*	`./` – positional argument specifying to recursively search the current working directory for appropriate files for `multiqc` to act on


**Input data:**

* the current working directory is specified, and `multiqc` recursively searches and grabs all appropriate files it can summarize

**Output data:**

* project_multiqc_output/project_multiqc_report.html (multiqc output html summary)
* project_multiqc_output/project_multiqc_data.zip (zipped directory containing multiqc output data)

<br>

---

## 11. Generate reference genome annotation information

### 11a. GTF to BED conversion
A bed-formatted annotation file is needed for adding annotation information to the output from the differential methylation analysis. We utilize gtf files from [Ensembl](https://www.ensembl.org/) and convert them as in the following example:

```bash
# downloading mouse reference gtf for this example
curl -LO https://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz

gunzip Mus_musculus.GRCm38.101.gtf.gz

gtfToGenePred Mus_musculus.GRCm38.101.gtf Mus_musculus.GRCm38.101.genePred

genePredToBed Mus_musculus.GRCm38.101.genePred Mus_musculus.GRCm38.101.bed

# removing intermediate file
rm Mus_musculus.GRCm38.101.genePred
```

**Input data:**

* a reference gtf file ("Mus_musculus.GRCm38.101.gtf" in the above example)

**Output data:**

* the generated bed file ("Mus_musculus.GRCm38.101.bed" in the above example)

### 11b. Making a mapping file of genes to transcripts
Making a mapping file of gene names to transcript names, which we need to link functional annotations in a primary output table. We can generate this map from the gtf file like so:

```bash
awk ' $3 == "transcript" ' Mus_musculus.GRCm38.101.gtf | cut -f 9 | tr -s ";" "\t" | \
    cut -f 1,3 | tr -s " " "\t" | cut -f 2,4 | tr -d '"' \
    > Mus_musculus.GRCm38.101-gene-to-transcript-map.tsv
```

**Input data:**

* a reference gtf file ("Mus_musculus.GRCm38.101.gtf" in the above example)

**Output data:**

* the generated file with gene IDs in the first column and transcript IDs in the second ("Mus_musculus.GRCm38.101-gene-to-transcript-map.tsv" in the above example)


### 11c. Making a table of all gene annotations
This is being extricated from the workflow and handled separately. In progress.


## 12. Differential methylation analysis

Example data for the R code below can be downloaded and unpacked with the following:

```bash
curl -L -o subset-test-results.tar https://figshare.com/ndownloader/files/34780726
tar -xvf subset-test-results.tar && rm subset-test-results.tar
```

The remainder of this section is performed in R. 

```R
library(tidyverse)
library(methylKit)

## reading in data
file.list <-list("subset-test-results/F-SRR12865062-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/F-SRR12865063-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/F-SRR12865064-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/G-SRR12865070-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/G-SRR12865071-sub_trimmed_bismark_bt2.bismark.cov.gz",
                 "subset-test-results/G-SRR12865072-sub_trimmed_bismark_bt2.bismark.cov.gz")

sample.list <- list("F-SRR12865062",
                    "F-SRR12865063",
                    "F-SRR12865064",
                    "G-SRR12865070",
                    "G-SRR12865071",
                    "G-SRR12865072")

# reading into memory
myobj <- methRead(location = file.list,
                  sample.id = sample.list,
                  assembly = "Mmus_GRCm39",
                  pipeline = "bismarkCoverage",
                  header = FALSE,
                  treatment = c(1,1,1,0,0,0),
                  mincov = 10)

# example of how to store as tables if memory requirements are too high
# myobj_storage <- methRead(location = file.list,
#                   sample.id = sample.list,
#                   assembly = "Mmus_GRCm39",
#                   dbtype = "tabix",
#                   pipeline = "bismarkCoverage",
#                   header = FALSE,
#                   treatment = c(1,1,1,0,0,0),
#                   dbdir = "methylkit-dbs/",
#                   mincov = 10)

## merging samples
meth <- unite(myobj)

## Finding diff methylated bases
myDiff <- calculateDiffMeth(meth, mc.cores = 4)

# get hyper methylated bases
myDiff25p.hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")
# making table for writing out
sig_hyper_out_tab <- getData(myDiff25p.hyper) %>% arrange(qvalue)

# get hypo methylated bases
myDiff25p.hypo <- getMethylDiff(myDiff, difference=25, qvalue = 0.01, type = "hypo")
# making table for writing out
sig_hypo_out_tab <- getData(myDiff25p.hypo) %>% arrange(qvalue)

# get all differentially methylated bases
myDiff25p <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01)
sig_all_out_tab <- getData(myDiff25p) %>% arrange(qvalue)

# writing out tables
write.table(sig_hyper_out_tab, "sig-hypermethylated-out.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_hypo_out_tab, "sig-hypomethylated-out.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_all_out_tab, "sig-all-methylated-out.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

## Annotating
library(genomation)
gene.obj <- readTranscriptFeatures("Mus_musculus.GRCm38.96.bed", up.flank = 1000, 
                                   down.flank = 1000, remove.unusual = TRUE, 
                                   unique.prom = TRUE)

diffAnn <- annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)

# making sig table with features 
sig_all_out_tab_with_features <- cbind(data.frame(myDiff25p), 
                                       getAssociationWithTSS(diffAnn), 
                                       as.data.frame(getMembers(diffAnn))) %>% .[,-c(8)]
write.table(sig_all_out_tab_with_features, "sig-all-methylated-out-with-features.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# making table of percent methylated
perc.meth <- percMethylation(meth, rowids = TRUE)
write.table(perc.meth, "percent-methylated.tsv", sep = "\t", 
            quote = FALSE, row.names = TRUE, col.names=NA)

# making sig table with features and functional annotations
## reading in annotation table appropriate for current organism
    ## when we have the final location for this information, this will
    ## need to be updated to pull a table with the links, rather than
    ## the link being hard-coded here in this example
functional_annots_tab <- 
    read.table("https://figshare.com/ndownloader/files/35939642", sep = "\t", 
               quote = "", header = TRUE)

# reading in gene to transcript mapping file
gene_transcript_map <- 
    read.table("Mus_musculus.GRCm38.101-gene-to-transcript-map.tsv", sep = "\t", 
               col.names = c("gene_ID", "feature.name"))

# for each transcript ID in the sig_all_out_tab_with_features table, getting 
# its corresponding gene ID and adding that to the table
sig_all_out_tab_with_features_and_gene_IDs <- 
    left_join(sig_all_out_tab_with_features, gene_transcript_map)

# now adding full annotations
sig_all_out_tab_with_features_and_annots <- 
    left_join(sig_all_out_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

# and writing out
write.table(sig_all_out_tab_with_features_and_annots, 
            "sig-all-methylated-out-with-features-and-annots.tsv", 
            sep = "\t", quote = FALSE, col.names=NA)
```

**Input data:**
* \*.bismark.cov.gz - gzip-compressed bedGraph-formatted files generated in Step 6 above
* Mus_musculus.GRCm38.101-gene-to-transcript-map.tsv - gene-to-transcript mapping file generated in Step 11b above

**Output data:**
* sig-hypermethylated-out.tsv - cytosines with significantly elevated methylation levels in treatment as compared to control
* sig-hypomethylated-out.tsv - cytosines with significantly reduced methylation levels in treatment as compared to control
* sig-all-methylated-out.tsv - all significantly differentially methylated cytosines
* sig-all-methylated-out-with-features.tsv - all significantly differentially methylated cytosines with features (promotor, exon, intron)
* sig-all-methylated-out-with-features-and-annots.tsv - all significantly differentially methylated cytosines with gene IDs, features (promotor, exon, intron), and functional annotations
* percent-methylated.tsv - table of methylation levels across all cytosines and samples
* sig-diff-meth-Cs-by-region.pdf - pie chart with percentages of differentially methylated cytosines

\* all of these files, except "percent-methylated.tsv", will be prefixed with contrasted groups, e.g. Group_1_vs_Group_2-\*

---
---
