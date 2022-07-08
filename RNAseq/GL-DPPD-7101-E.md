# GeneLab bioinformatics processing pipeline for Illumina RNA-sequencing data

> **This page holds an overview and instructions for how GeneLab processes RNAseq datasets. Exact processing commands and GL-DPPD-7101 version used for specific datasets are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and processed data output files are provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** April 26, 2022  
**Revision:** E  
**Document Number:** GL-DPPD-7101-E  

**Submitted by:**  
Jonathan Oribello (GeneLab Data Processing Team)

**Approved by:**  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Jonathan Galazka (GeneLab Project Scientist)

---

## Updates from previous version

[Software used](#software-used) now specifies exact version numbers.

The following multiQC compilation steps now force interactive plots: 
- [step 1b](#1b-compile-raw-data-qc), [step 2c](#2c-compile-trimmed-data-qc), [step 4b](#4b-compile-alignment-logs), [step 6b](#6b-compile-strandedness-reports).

Updated [Ensembl Reference Files](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv) now used:
- Animals: Ensembl release 101
- Plants: Ensembl plants release 48
- Bacteria: Ensembl bacteria release 48

STAR Gene Counts now generated in [step 4a](#4a-align-reads-to-reference-genome-with-star)

- These counts are tabulated in new [step 4c](#4c-tablulate-star-counts-in-r)

Aligned reads are now subsequently sorted with Samtools in new [step 4d](#4d-sort-aligned-reads) 
- Sorted aligned reads are then indexed with Samtools in new [step 4e](#4e-index-sorted-aligned-reads)

Additional RSeQC analyses are performed on genome aligned reads as follows:

- GeneBody coverage is evaluated and reports are compiled with multiQC in [step 6c](#6c-evaluate-genebody-coverage) and [step 6d](#6d-compile-genebody-coverage-reports), respectively
- For paired end datasets, inner distance is determined and reports are compiled with multiQC in [step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only) and [step 6f](#6f-compile-inner-distance-reports), respectively
- Read distribution is assessed and reports are compiled with multiQC in [step 6g](#6g-assess-read-distribution) and [step 6h](#6h-compile-read-distribution-reports), respectively

RSEM now quantitates all reads >= 20bp
- The `--seed-length 20` option was added in [step 8a](#8a-count-aligned-reads-with-rsem)

RSEM Count results are additionally summarized as follows:

- MultiQC is used to compile RSEM count reports in [step 8b](#8b-compile-rsem-count-logs)
- The total number of genes expressed per sample are calculated in new [step 8c](#8c-calculate-total-number-of-genes-expressed-per-sample-in-r)

The DESeq2 Normalization and DGE step for datasets with ERCC spike-in, [step 9a](#9a-for-datasets-with-ercc-spike-in), was modified as follows:

- Fixed bug where `ERCCnorm_contrasts.csv` was always the same as the non-ERCC contrasts.csv
  > Note: In most cases, these files are the same. They will only differ when, for the ERCC-based analysis, removal of samples with no detectable Group B ERCC spike-in results in the a complete removal of a group.

The DESeq2 Normalization and DGE step for both datasets with ERCC spike-in, [step 9a](#9a-for-datasets-with-ercc-spike-in), and without, [step 9b](#9b-for-datasets-without-ercc-spike-in) was modified as follows:

- Input file regex modified to address bug that occured when certain sample IDs were substrings of other sample IDs (e.g. Sample1, Sample13)
- Output file: `Unnormalized_Counts.csv` renamed to `RSEM_Unnormalized_Counts.csv` for clarity

ERCC Analysis is performed as detailed in [step 10](#10-evaluate-ercc-spike-in-data):

- ERCC Counts are plotted and quantified [step 10a](#10a-evaluate-ercc-count-data-in-python)
- DESeq2 differential gene expression in Mix 1 versus Mix 2 groups is performed on ERCC counts [step 10b](#10b-perform-deseq2-analysis-of-ercc-counts-in-r)
- DESeq2 results are analyzed for the expected Mix 1 vs Mix 2 ERCC groups A, B, C, and D ratios [step 10c](#10c-analyze-ercc-deseq2-results-in-python)

---

# Table of contents  

- [GeneLab bioinformatics processing pipeline for Illumina RNA-sequencing data](#genelab-bioinformatics-processing-pipeline-for-illumina-rna-sequencing-data)
  - [Updates from previous version](#updates-from-previous-version)
- [Table of contents](#table-of-contents)
- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Raw Data QC](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [2. Trim/Filter Raw Data and Trimmed Data QC](#2-trimfilter-raw-data-and-trimmed-data-qc)
    - [2a. Trim/Filter Raw Data](#2a-trimfilter-raw-data)
    - [2b. Trimmed Data QC](#2b-trimmed-data-qc)
    - [2c. Compile Trimmed Data QC](#2c-compile-trimmed-data-qc)
  - [3. Build STAR Reference](#3-build-star-reference)
  - [4. Align Reads to Reference Genome then Sort and Index](#4-align-reads-to-reference-genome-then-sort-and-index)
    - [4a. Align Reads to Reference Genome with STAR](#4a-align-reads-to-reference-genome-with-star)
    - [4b. Compile Alignment Logs](#4b-compile-alignment-logs)
    - [4c. Tablulate STAR Counts in R](#4c-tablulate-star-counts-in-r)
    - [4d. Sort Aligned Reads](#4d-sort-aligned-reads)
    - [4e. Index Sorted Aligned Reads](#4e-index-sorted-aligned-reads)
  - [5. Create Reference BED File](#5-create-reference-bed-file)
    - [5a. Convert GTF to genePred File](#5a-convert-gtf-to-genepred-file)
    - [5b. Convert genePred to BED File](#5b-convert-genepred-to-bed-file)
  - [6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC](#6-assess-strandedness-genebody-coverage-inner-distance-and-read-distribution-with-rseqc)
    - [6a. Determine Read Strandedness](#6a-determine-read-strandedness)
    - [6b. Compile Strandedness Reports](#6b-compile-strandedness-reports)
    - [6c. Evaluate GeneBody Coverage](#6c-evaluate-genebody-coverage)
    - [6d. Compile GeneBody Coverage Reports](#6d-compile-genebody-coverage-reports)
    - [6e. Determine Inner Distance (For Paired End Datasets ONLY)](#6e-determine-inner-distance-for-paired-end-datasets-only)
    - [6f. Compile Inner Distance Reports](#6f-compile-inner-distance-reports)
    - [6g. Assess Read Distribution](#6g-assess-read-distribution)
    - [6h. Compile Read Distribution Reports](#6h-compile-read-distribution-reports)
  - [7. Build RSEM Reference](#7-build-rsem-reference)
  - [8. Quantitate Aligned Reads](#8-quantitate-aligned-reads)
    - [8a. Count Aligned Reads with RSEM](#8a-count-aligned-reads-with-rsem)
    - [8b. Compile RSEM Count Logs](#8b-compile-rsem-count-logs)
    - [8c. Calculate Total Number of Genes Expressed Per Sample in R](#8c-calculate-total-number-of-genes-expressed-per-sample-in-r)
  - [9. Create Runsheet](#9-create-runsheet)
    - [9a. Download ISA Archive](#9a-download-isa-archive)
    - [9b. Generate runsheet from ISA Archive](#9b-generate-runsheet-from-isa-archive)
  - [10. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R](#10-normalize-read-counts-perform-differential-gene-expression-analysis-and-add-gene-annotations-in-r)
    - [10a. Using default 'median of ratios method' Deseq2 normalization](#10a-using-default-median-of-ratios-method-deseq2-normalization)
    - [10b. Using Deseq2 ERCC group B based normalization](#10b-using-deseq2-ercc-group-b-based-normalization)
  - [11. Evaluate ERCC Spike-In Data](#11-evaluate-ercc-spike-in-data)
    - [11a. Evaluate ERCC Count Data in Python](#11a-evaluate-ercc-count-data-in-python)
    - [11b. Perform DESeq2 Analysis of ERCC Counts in R](#11b-perform-deseq2-analysis-of-ercc-counts-in-r)
    - [11c. Analyze ERCC DESeq2 Results in Python](#11c-analyze-ercc-deseq2-results-in-python)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.11.9|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.12|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|3.7|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|0.6.7|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|2.7.10a|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|RSEM|1.3.1|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Samtools|1.15|[http://www.htslib.org/](http://www.htslib.org/)|
|gtfToGenePred|377|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed|377|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|infer_experiment|4.0.0|[http://rseqc.sourceforge.net/#infer-experiment-py](http://rseqc.sourceforge.net/#infer-experiment-py)|
|geneBody_coverage|4.0.0|[http://rseqc.sourceforge.net/#genebody-coverage-py](http://rseqc.sourceforge.net/#genebody-coverage-py)|
|inner_distance|4.0.0|[http://rseqc.sourceforge.net/#inner-distance-py](http://rseqc.sourceforge.net/#inner-distance-py)|
|read_distribution|4.0.0|[http://rseqc.sourceforge.net/#read-distribution-py](http://rseqc.sourceforge.net/#read-distribution-py)|
|R|4.1.2|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.14.0|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|1.34|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|1.22|[https://bioconductor.org/packages/release/bioc/html/tximport.html](https://bioconductor.org/packages/release/bioc/html/tximport.html)|
|tidyverse|1.3.1|[https://www.tidyverse.org](https://www.tidyverse.org)|
|Risa|1.36|[https://www.bioconductor.org/packages/release/bioc/html/Risa.html](https://www.bioconductor.org/packages/release/bioc/html/Risa.html)|

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are provided in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory.
> 
> All output files marked with a \# are published for each RNAseq processed dataset in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 

---

## 1. Raw Data QC

<br>

### 1a. Raw Data QC  

```bash
fastqc -o /path/to/raw_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

### 1b. Compile Raw Data QC  

```bash
multiqc --interactive -n raw_multiqc -o /path/to/raw_multiqc/output/directory /path/to/directory/containing/raw_fastqc/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

- raw_multiqc.html\# (multiqc report)
- /raw_multiqc_data\# (directory containing multiqc data)

<br>

---

## 2. Trim/Filter Raw Data and Trimmed Data QC

<br>

### 2a. Trim/Filter Raw Data  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores NumberOfThreads \
  --phred33 \
  --illumina \ # if adapters are not illumina, replace with adapters used
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \ # only for PE studies, remove this paramater if raw data are SE
  sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz
# if SE, replace the last line with only the forward reads (R1) of each sample

```

**Parameter Definitions:**

- `--gzip` – compress the output files with `gzip`
- `--path_to_cutadapt` - specify path to cutadapt software if it is not in your `$PATH`
- `--cores` - specify the number of threads available on the server node to perform trimming
- `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
- `--illumina` - defines the adapter sequence to be trimmed as the first 13bp of the Illumina universal adapter `AGATCGGAAGAGC`
- `--output_dir` - the output directory to store results
- `--paired` - indicates paired-end reads - both reads, forward (R1) and reverse (R2) must pass length threshold or else both reads are removed
- `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- *fastq.gz\# (trimmed reads)
- *trimming_report.txt\# (trimming report)

<br>

### 2b. Trimmed Data QC  

```bash
fastqc -o /path/to/trimmed_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**

- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

### 2c. Compile Trimmed Data QC  

```bash
multiqc --interactive -n trimmed_multiqc -o /path/to/trimmed_multiqc/output/directory /path/to/directory/containing/trimmed_fastqc/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/trimmed_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 2b](#2b-trimmed-data-qc))

**Output Data:**

- trimmed_multiqc.html\# (multiqc report)
- /trimmed_multiqc_data\# (directory containing multiqc data)

<br>

---

## 3. Build STAR Reference  

```bash
STAR --runThreadN NumberOfThreads \
  --runMode genomeGenerate \
  --limitGenomeGenerateRAM 55000000000 \
  --genomeSAindexNbases 14 \
  --genomeDir /path/to/STAR/genome/directory \
  --genomeFastaFiles /path/to/genome/fasta/file \
  --sjdbGTFfile /path/to/annotation/gtf/file \
  --sjdbOverhang ReadLength-1

```

**Parameter Definitions:**

- `--runThreadN` – number of threads available on server node to create STAR reference
- `--runMode` - instructs STAR to run genome indices generation job
- `--limitGenomeGenerateRAM` - maximum RAM available (in bytes) to generate STAR reference, at least 35GB are needed for mouse and the example above shows 55GB
- `--genomeSAindexNbases` - length (in bases) of the SA pre-indexing string, usually between 10 and 15. Longer strings require more memory but allow for faster searches. This value should be scaled down for smaller genomes (like bacteria) to min(14, log2(GenomeLength)/2 - 1). For example, for a 1 megaBase genome this value would be 9.
- `--genomeDir` - specifies the path to the directory where the STAR reference will be stored. At least 100GB of available disk space is required for mammalian genomes.
- `--genomeFastaFiles` - specifies one or more fasta file(s) containing the genome reference sequences
- `--sjdbGTFfile` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `--sjdbOverhang` - indicates the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. The length should be one less than the maximum length of the reads.

**Input Data:**

- *.fasta ([genome sequence](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))
- *.gtf ([genome annotation](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))

**Output Data:**

STAR genome reference, which consists of the following files:

- chrLength.txt
- chrNameLength.txt
- chrName.txt
- chrStart.txt
- exonGeTrInfo.tab
- exonInfo.tab
- geneInfo.tab
- Genome
- genomeParameters.txt
- SA
- SAindex
- sjdbInfo.txt
- sjdbList.fromGTF.out.tab
- sjdbList.out.tab
- transcriptInfo.tab

<br>

---

## 4. Align Reads to Reference Genome then Sort and Index

<br>

### 4a. Align Reads to Reference Genome with STAR

```bash
STAR --twopassMode Basic \
 --limitBAMsortRAM 65000000000 \
 --genomeDir /path/to/STAR/genome/directory \
 --outSAMunmapped Within \
 --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD MC \
 --outFilterMultimapNmax 20 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \ # for PE only
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --sjdbScore 1 \
 --readFilesCommand zcat \
 --runThreadN NumberOfThreads \
 --outSAMtype BAM SortedByCoordinate \
 --quantMode TranscriptomeSAM GeneCounts \
 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
 --outFileNamePrefix /path/to/STAR/output/directory/<sample_id> \
 --readFilesIn /path/to/trimmed_forward_reads \
 /path/to/trimmed_reverse_reads # only needed for PE studies

```

**Parameter Definitions:**

- `--twopassMode` – specifies 2-pass mapping mode; the `Basic` option instructs STAR to perform the 1st pass mapping, then automatically extract junctions, insert them into the genome index, and re-map all reads in the 2nd mapping pass
- `--limitBAMsortRAM` - maximum RAM available (in bytes) to sort the bam files, the example above indicates 65GB
- `--genomeDir` - specifies the path to the directory where the STAR reference is stored
- `--outSAMunmapped` - specifies ouput of unmapped reads in the sam format; the `Within` option instructs STAR to output the unmapped reads within the main sam file
- `--outFilterType` - specifies the type of filtering; the `BySJout` option instructs STAR to keep only those reads that contain junctions that passed filtering in the SJ.out.tab output file
- `--outSAMattributes` - list of desired sam attributes in the order desired for the output sam file; sam attribute descriptions can be found [here](https://samtools.github.io/hts-specs/SAMtags.pdf)
- `--outFilterMultimapNmax` – specifies the maximum number of loci the read is allowed to map to; all alignments will be output only if the read maps to no more loci than this value
- `--outFilterMismatchNmax` - maximum number of mismatches allowed to be included in the alignment output
- `--outFilterMismatchNoverReadLmax` - ratio of mismatches to read length allowed to be included in the alignment output; the `0.04` value indicates that up to 4 mismatches are allowed per 100 bases
- `--alignIntronMin` - minimum intron size; a genomic gap is considered an intron if its length is equal to or greater than this value, otherwise it is considered a deletion
- `--alignIntronMax` - maximum intron size
- `--alignMatesGapMax` - maximum genomic distance (in bases) between two mates of paired-end reads; this option should be removed for single-end reads
- `--alignSJoverhangMin` - minimum overhang (i.e. block size) for unannotated spliced alignments
- `--alignSJDBoverhangMin` - minimum overhang (i.e. block size) for annotated spliced alignments
- `--sjdbScore` - additional alignment score for alignments that cross database junctions
- `--readFilesCommand` - specifies command needed to interpret input files; the `zcat` option indicates input files are compressed with gzip and zcat will be used to uncompress the gzipped input files
- `--runThreadN` - indicates the number of threads to be used for STAR alignment and should be set to the number of available cores on the server node
- `--outSAMtype` - specifies desired output format; the `BAM SortedByCoordinate` options specify that the output file will be sorted by coordinate and be in the bam format
- `--quantMode` - specifies the type(s) of quantification desired; the `TranscriptomeSAM` option instructs STAR to output a separate sam/bam file containing alignments to the transcriptome and the `GeneCounts` option instructs STAR to output a tab delimited file containing the number of reads per gene
- `--outSAMheaderHD` - indicates a header line for the sam/bam file
- `--outFileNamePrefix` - specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id
- `--readFilesIn` - path to input read 1 (forward read) and read 2 (reverse read); for paired-end reads, read 1 and read 2 should be separated by a space; for single-end reads only read 1 should be indicated

**Input Data:**

- STAR genome reference (output from [Step 3](#3-build-star-reference))
- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam\# (sorted mapping to transcriptome)
- *Log.final.out\# (log file containing alignment info/stats such as reads mapped, etc)
- *ReadsPerGene.out.tab (tab deliminated file containing STAR read counts per gene with 4 columns that correspond to different strandedness options: column 1 = gene ID, column 2 = counts for unstranded RNAseq, column 3 = counts for 1st read strand aligned with RNA, column 4 = counts for 2nd read strand aligned with RNA)
- *Log.out
- *Log.progress.out
- *SJ.out.tab\# (high confidence collapsed splice junctions in tab-delimited format)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

<br>

### 4b. Compile Alignment Logs

```bash
multiqc --interactive -n align_multiqc -o /path/to/aligned_multiqc/output/directory /path/to/*Log.final.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*Log.final.out/files` – the directory holding the *Log.final.out output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc., output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- align_multiqc.html\# (multiqc report)
- /align_multiqc_data\# (directory containing multiqc data)

<br>

### 4c. Tablulate STAR Counts in R

```R
print("Make STAR counts table")
print("")

work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain samples.txt file
align_dir="/path/to/directory/containing/STAR/counts/files"

setwd(file.path(work_dir))

### Pull in sample names where the "samples.txt" file is a single column list of sample names ###
study <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import Data
ff <- list.files(file.path(align_dir), pattern = "ReadsPerGene.out.tab", recursive=TRUE, full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the ISA samples
ff <- ff[sapply(rownames(study), function(x)grep(paste0(x,'_ReadsPerGene.out.tab$'), ff, value=FALSE))]

# Remove the first 4 lines
counts.files <- lapply( ff, read.table, skip = 4 )

# Get counts aligned to either strand for unstranded data by selecting col 2, to the first (forward) strand by selecting col 3 or to the second (reverse) strand by selecting col 4
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 3 ] ) )

# Add column and row names
colnames(counts) <- rownames(study)
row.names(counts) <- counts.files[[1]]$V1


##### Export unnormalized counts table
setwd(file.path(align_dir))
write.csv(counts,file='STAR_Unnormalized_Counts.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (A newline delimited list of sample IDs)
- *ReadsPerGene.out.tab (STAR counts per gene, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- STAR_Unnormalized_Counts.csv\# (Table containing raw STAR counts for each sample)

<br>

### 4d. Sort Aligned Reads

```bash
samtools sort -m 3G \
	--threads NumberOfThreads \
	-o /path/to/*Aligned.sortedByCoord_sorted.out.bam \
  /path/to/*Aligned.sortedByCoord.out.bam
```

**Parameter Definitions:**

- `-m` - memory available per thread, `3G` indicates 3 gigabytes, this can be changed based on user resources
- `--threads` - number of threads available on server node to sort genome alignment files
- `/path/to/*Aligned.sortedByCoord.out.bam` – path to the *Aligned.sortedByCoord.out.bam output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome file, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- *Aligned.sortedByCoord_sorted.out.bam\# (samtools sorted genome aligned bam file)

<br>

### 4e. Index Sorted Aligned Reads

```bash
samtools index -@ NumberOfThreads /path/to/*Aligned.sortedByCoord_sorted.out.bam
```

**Parameter Definitions:**

- `-@` - number of threads available on server node to index the sorted alignment files
- `/path/to/*Aligned.sortedByCoord_sorted.out.bam` – the path to the sorted *Aligned.sortedByCoord_sorted.out.bam output files from the [step 4d](#4d-sort-aligned-reads), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, ourput from [Step 4d](#4d-sort-aligned-reads))

**Output Data:**

- *Aligned.sortedByCoord_sorted.out.bam.bai\# (index of sorted mapping to genome file)

<br>

---

## 5. Create Reference BED File

<br>

### 5a. Convert GTF to genePred File  

```bash
gtfToGenePred /path/to/annotation/gtf/file \
  /path/to/output/genePred/file

```

**Parameter Definitions:**

- `/path/to/annotation/gtf/file` – specifies the file(s) containing annotated reference transcripts in the standard gtf format, provided as a positional argument
- `/path/to/output/genePred/file` – specifies the location and name of the output genePred file(s), provided as a positional argument

**Input Data:**

- *.gtf ([genome annotation](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))

**Output Data:**

- *.genePred (genome annotation in genePred format)

<br>

### 5b. Convert genePred to BED File  

```bash
genePredToBed /path/to/annotation/genePred/file \
  /path/to/output/BED/file
```

**Parameter Definitions:**

- `/path/to/annotation/genePred/file` – specifies the file(s) containing annotated reference transcripts in the genePred format, provided as a positional argument
- `/path/to/output/BED/file` – specifies the location and name of the output BED file(s), provided as a positional argument

**Input Data:**

- *.genePred (genome annotation in genePred format, output from [Step 5a](#5a-convert-gtf-to-genepred-file))

**Output Data:**

- *.bed (genome annotation in BED format)

<br>

---

## 6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC

<br>

### 6a. Determine Read Strandedness

```bash
infer_experiment.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -s 15000000 > /path/to/*infer_expt.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-s` - specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `>` - redirects standard output to specified file
- `/path/to/*infer_expt.out` - specifies the location and name of the file containing the infer_experiment standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *infer_expt.out (file containing the infer_experiment standard output)

<br>

### 6b. Compile Strandedness Reports

```bash
multiqc --interactive -n infer_exp_multiqc -o /path/to/infer_exp_multiqc/output/directory /path/to/*infer_expt.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*infer_expt.out/files` – the directory holding the *infer_expt.out output files from the [read strandedness step](#6a-determine-read-strandedness), provided as a positional argument

**Input Data:**

- *infer_expt.out (file containing the infer_experiment standard output, output from [Step 6a](#6a-determine-read-strandedness))

**Output Data:**

- infer_exp_multiqc.html\# (multiqc report)
- /infer_exp_multiqc_data\# (directory containing multiqc data)

<br>

### 6c. Evaluate GeneBody Coverage

```bash
geneBody_coverage.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -o /path/to/geneBody_coverage/output/directory/<sample_id>
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-o` - specifies the path to the output directory
- `/path/to/geneBody_coverage/output/directory/<sample_id>` - specifies the location and name of the directory containing the geneBody_coverage output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.geneBodyCoverage.curves.pdf (genebody coverage line plot)
- *.geneBodyCoverage.r (R script that generates the genebody coverage line plot)
- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values used to generate the line plot)

<br>

### 6d. Compile GeneBody Coverage Reports

```bash
multiqc --interactive -n genebody_cov_multiqc -o /path/to/geneBody_coverage_multiqc/output/directory /path/to/geneBody_coverage/output/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/geneBody_coverage/output/files` – the directory holding the geneBody_coverage output files from [step 6c](#6c-evaluate-genebody-coverage), provided as a positional argument

**Input Data:**

- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values, output from [Step 6c](#6c-evaluate-genebody-coverage))

**Output Data:**

- geneBody_cov_multiqc.html\# (multiqc report)
- /geneBody_cov_multiqc_data\# (directory containing multiqc data)

<br>

### 6e. Determine Inner Distance (For Paired End Datasets ONLY)

```bash
inner_distance.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -k 15000000 \
 -l -150 \
 -u 350 \
 -o  /path/to/inner_distance/output/directory
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-k` - specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `-l` - specifies the lower bound of inner distance (bp).
- `-u` - specifies the upper bound of inner distance (bp)
- `/path/to/inner_distance/output/directory` - specifies the location and name of the directory containing the inner_distance output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.inner_distance.txt (log of read-wise inner distance results)
- *.inner_distance_freq.txt (tab delimited table of inner distances mapped to number of reads with that distance)
- *.inner_distance_plot.pdf (histogram plot of inner distance distribution)
- *.inner_distance_plot.r (R script that generates the histogram plot)

<br>

### 6f. Compile Inner Distance Reports

```bash
multiqc --interactive -n inner_dist_multiqc /path/to/inner_dist_multiqc/output/directory /path/to/inner_dist/output/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/inner_dist/output/files` – the directory holding the inner_distance output files from [Step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only), provided as a positional argument

**Input Data:**

- *.inner_distance_freq.txt (tab delimited table of inner distances from [step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only))

**Output Data:**

- inner_distance_multiqc.html\# (multiqc report)
- /inner_distance_multiqc_data\# (directory containing multiqc data)

<br>

### 6g. Assess Read Distribution

```bash
read_distribution.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam > /path/to/*read_dist.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `>` - redirects standard output to specified file
- `/path/to/*read_dist.out` - specifies the location and name of the file containing the read_distribution standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *read_dist.out (file containing the read distribution standard output)

<br>

### 6h. Compile Read Distribution Reports

```bash
multiqc --interactive -n read_dist_multiqc -o /path/to/read_dist_multiqc/output/directory /path/to/*read_dist.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*read_dist.out/files` – the directory holding the *read_dist.out output files from [Step 6g](#6g-assess-read-distribution) provided as a positional argument

**Input Data:**

- *read_dist.out (files containing the read_distributation standard output, output from [Step 6g](#6g-assess-read-distribution))

**Output Data:**

- read_dist_multiqc.html\# (multiqc report)
- /read_dist_multiqc_data\# (directory containing multiqc data)

<br>

---

## 7. Build RSEM Reference

```bash
rsem-prepare-reference --gtf /path/to/annotation/gtf/file \
 /path/to/genome/fasta/file \
 /path/to/RSEM/genome/directory/RSEM_ref_prefix

```

**Parameter Definitions:**

- `--gtf` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `/path/to/genome/fasta/file` – specifies one or more fasta file(s) containing the genome reference sequences, provided as a positional argument
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference will be stored and the prefix desired for the RSEM reference files, provided as a positional argument

**Input Data:**

- *.fasta ([genome sequence](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))
- *.gtf ([genome annotation](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))

**Output Data:**

RSEM genome reference, which consists of the following files:

- RSEM_ref_prefix.chrlist
- RSEM_ref_prefix.grp
- RSEM_ref_prefix.idx.fa
- RSEM_ref_prefix.n2g.idx.fa
- RSEM_ref_prefix.seq
- RSEM_ref_prefix.ti
- RSEM_ref_prefix.transcripts.fa

<br>

---

## 8. Quantitate Aligned Reads

<br>

### 8a. Count Aligned Reads with RSEM

```bash
rsem-calculate-expression --num-threads NumberOfThreads \
 --alignments \
 --bam \
 --paired-end \
 --seed 12345 \
 --seed-length 20 \
 --estimate-rspd \
 --no-bam-output \
 --strandedness reverse|forward|none \
 /path/to/*Aligned.toTranscriptome.out.bam \
 /path/to/RSEM/genome/directory/RSEM_ref_prefix \
 /path/to/RSEM/counts/output/directory/<sample_id>
```

**Parameter Definitions:**

- `--num-threads` – specifies the number of threads to use
- `--alignments` - indicates that the input file contains alignments in sam, bam, or cram format
- `--bam` - specifies that the input alignments are in bam format
- `--paired-end` - indicates that the input reads are paired-end reads; this option should be removed if the input reads are single-end
- `--seed` - the seed for the random number generators used in calculating posterior mean estimates and credibility intervals; must be a non-negative 32-bit integer
- `--seed-length 20` - instructs RSEM to ignore any aligned read if it or its mates' (for paired-end reads) length is less than 20bp
- `--estimate-rspd` - instructs RSEM to estimate the read start position distribution (rspd) from the data
- `--no-bam-output` - instructs RSEM not to output any bam file
- `--strandedness` - defines the strandedness of the RNAseq reads; the `reverse` option is used if read strandedness (output from [step 6](#6a-determine-read-strandedness)) is antisense, `forward` is used with sense strandedness, and `none` is used if strandedness is half sense half antisense
- `/path/to/*Aligned.toTranscriptome.out.bam` - specifies path to input bam files, provided as a positional argument
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference is stored and its prefix, provided as a positional argument
- `/path/to/RSEM/counts/output/directory` – specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id

**Input Data:**

- RSEM genome reference (output from [Step 7](#7-build-rsem-reference))
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- *genes.results\# (counts per gene)
- *isoforms.results\# (counts per isoform)
- *stat (directory containing the following stats files)
  - *cnt
  - *model
  - *theta

<br>

### 8b. Compile RSEM Count Logs

```bash
multiqc --interactive -n RSEM_count_multiqc -o /path/to/RSEM_count_multiqc/output/directory /path/to/*stat/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*stat/files` – the directories holding the *stat output files from the [RSEM Counts step](#8a-count-aligned-reads-with-rsem), provided as a positional argument

**Input Data:**

- *stat (directory containing the following stats files, output from [Step 8a](#8a-count-aligned-reads-with-rsem))
  - *cnt
  - *model
  - *theta

**Output Data:**

- RSEM_count_multiqc.html\# (multiqc report)
- /RSEM_count_multiqc_data\# (directory containing multiqc data)

<br>

### 8c. Calculate Total Number of Genes Expressed Per Sample in R

```R
library(tximport)
library(tidyverse)

work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain samples.txt file
counts_dir="/path/to/directory/containing/RSEM/counts/files"

setwd(file.path(work_dir))

### Pull in sample names where the "samples.txt" file is a single column list of sample names ###
samples <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import RSEM Gene Count Data
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)

### reorder the genes.results files to match the ordering of the samples in the metadata file
files <- files[sapply(rownames(samples), function(x)grep(paste0(x,'.genes.results$'), files, value=FALSE))]

names(files) <- rownames(samples)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

##### Count the number of genes with non-zero counts for each sample 
rawCounts <- txi.rsem$counts
NumNonZeroGenes <- (as.matrix(colSums(rawCounts > 0), row.names = 1))
colnames(NumNonZeroGenes) <- c("Number of genes with non-zero counts")

##### Export the number of genes with non-zero counts for each sample
setwd(file.path(counts_dir))
write.csv(NumNonZeroGenes,file='NumNonZeroGenes.csv')

## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (A newline delimited list of sample IDs)
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))

**Output Data:**

- NumNonZeroGenes.csv (A samplewise table of the number of genes expressed)

<br>

---

## 9. Create Runsheet

The runsheet is a csv file that includes sample metadata required for running the R scripts in [Step 10](TBA).  For GLDS datasets, this can be generated from the GLDS API using the dp_tools library cli tool from [github](https://github.com/J-81/dp_tools/releases/tag/1.0.7rc2). For non-GLDS datasets, this runsheet can be manually created, an example runsheet is located TBA[here]()

### 9a. Download ISA Archive

Note: An alternative to using this script is to download the ISA archive from the GLDS repository directly.

```bash
dpt-get-isa-archive --accession <gldsAccession> --alternate-url
```
**Parameter Definitions:**

- `--accession` - Required argument flag for supplying GLDS accession ID
- `gldsAccession` - GLDS accession ID to generate a runsheet
- `--alternate-url` - Use alternative and more stable GLDS api uri

**Input:**

- No input files for this step

**Output:**

- *.zip - (ISA Archive zip file for GLDS accession ID queried)

### 9b. Generate runsheet from ISA Archive

```bash
dpt-isa-to-runsheet --accession <gldsAccession> --config-type bulkRNASeq --isa-archive *.zip
```

**Parameter Definitions:**

- `--accession` - Required argument flag for supplying GLDS accession ID
- `gldsAccession` - GLDS accession ID to generate a runsheet
- `--config-type` - Required argument flag for specifying assay type to extract appropriatte metadata for assay type processing
- `--isa-archive` - ISA Archive zip file to extract metadata to generate the runsheet

**Input:**

- *.zip - (ISA Archive zip file for GLDS accession ID queried)

**Output:**

- *_bulkRNASeq_v1_runsheet.csv - (A csv file describing experimental groups and metadata required for scripts in [Step 10](TBA). If generated by using the dp_tools script, this also includes additional information used for the automated Nextflow implementation of this workflow)

---

## 10. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R

Code used in this step in located with the workflow codebase bin [directory](workflow_code/Nextflow_RCP/bin/dge_annotation_R_scripts)
The gene annotations are sourced from annotation database tables generated in the GeneLab [Annotation Database Table Generation workflow](../Annotation_Database_Table_Generation).
<br>

### 10a. Using default 'median of ratios method' Deseq2 normalization 

```R
./dge_annotation_R_scripts/dge_annotation_workflow.R \
    --runsheet_path <path/to/runsheet.csv> \
    --input_gene_results_dir <path/to/gene_results_directory> \
    --primary_keytype <TAIR|ENSEMBL> \
    --organisms_csv <organisms_csv> \
    --organism <organism_non_sci> \
    --normalization 'default' \
    --normalized_counts_output_prefix "norm_counts_output/" \
    --dge_output_prefix "dge_output/" \
    --extended_table_output_prefix "dge_output/"\
    --extended_table_output_suffix ".csv"
```

**Parameter Definitions:**

- `--runsheet_path` - Flag to specify the runsheet location
- `--input_gene_results_dir` - Flag to specify the gene.results files containing directory location
- `<path/to/gene_results_directory>` - Directory containing all *.gene.results files
- `--primary_keytype` - Annotation key used for joining additional gene annotations, currently supports 'TAIR' for Arobidopsis Thaliana and 'ENSEMBL' for the following organisms: 'Mus Musculus, Drosophila melanogaster, Homo sapiens'
- `--organisms_csv` - Flag to specify the location of the organisms.csv file that includes mappings to the organism specific gene annotation database tables URIs
- `--organism` - Flag to denote the row of organisms.csv to use
- `--normalization` - Flag to choose the normalization method. Supports the following options: 'default' - use default Deseq2 median of ratios approach using all unfiltered non-ERCC gene counts, 'ERCC-groupB' - uses ERCC group B spike-in gene counts as controlGroup genes. See Deseq2 manual at this [section](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#control-features-for-estimating-size-factors) for additional details.
- `--normalized_counts_output_prefix` - Flag to specify the location of the normalized counts output files
- `--dge_output_prefix` - Flag to specify the location of the differential gene expression out files
- `--extended_table_output_prefix` - Flag to specify the extended table output files prefix (Note: these files are used internally by GeneLab for the visualization portal)
- `--extended_table_output_suffix` - Flag to specify the extended table output files suffix (Note: these files are used internally by GeneLab for the visualization portal)

**Input Data:**

- path/to/runsheet.csv (runsheet csv file generated in [Step 9](#9-create-runsheet))
- [organisms.csv](https://github.com/J-81/Nextflow_RCP/blob/dev_rc1.0.6/assets/organisms.csv) (csv file containing short name, species name, taxon ID, and annotation db object of model organisms hosted on GeneLab)
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))

**Output Data:**

- RSEM_Unnormalized_Counts.csv\# (table containing raw RSEM gene counts for each sample)
- Normalized_Counts.csv\# (table containing normalized gene counts for each sample)
- SampleTable.csv\# (table containing samples and their respective groups)
- visualization_output_table.csv (file used to generate GeneLab DGE visualizations)
- visualization_PCA_table.csv (file used to generate GeneLab PCA plots)
- differential_expression.csv\# (table containing normalized counts for each sample, group statistics, DESeq2 DGE results for each pairwise comparison, and gene annotations) 
- contrasts.csv\# (table containing all pairwise comparisons)

<br>

---

### 10b. Using Deseq2 ERCC group B based normalization 

```R
./dge_annotation_R_scripts/dge_annotation_workflow.R \
    --runsheet_path runsheet.csv \
    --input_gene_results_dir <path/to/gene_results_directory> \
    --primary_keytype <TAIR|ENSEMBL> \
    --organisms_csv <organisms_csv> \
    --organism <organism_non_sci> \
    --normalization 'ERCC-groupB' \
    --normalized_counts_output_prefix "norm_counts_output/ERCC_" \
    --dge_output_prefix "dge_output_ercc/ERCCnorm_" \
    --extended_table_output_prefix "dge_output_ercc/"\
    --extended_table_output_suffix "_ERCCnorm.csv" \
    --verbose
```

**Parameter Definitions:**

- `--runsheet_path` - Flag to specify the runsheet location
- `--input_gene_results_dir` - Flag to specify the gene.results files containing directory location
- `<path/to/gene_results_directory>` - Directory containing all *.gene.results files
- `--primary_keytype` - Annotation key used for joining additional gene annotations, currently supports 'TAIR' for Arobidopsis Thaliana and 'ENSEMBL' for the following organisms: 'Mus Musculus, Drosophila melanogaster, Homo sapiens'
- `--organisms_csv` - Flag to specify the location of the organisms.csv file that includes mappings to the organism specific gene annotation database tables URIs
- `--organism` - Flag to denote the row of organisms.csv to use
- `--normalization` - Flag to choose the normalization method. Supports the following options: 'default' - use default Deseq2 median of ratios approach using all unfiltered non-ERCC gene counts, 'ERCC-groupB' - uses ERCC group B spike-in gene counts as controlGroup genes. See Deseq2 manual at this [section](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#control-features-for-estimating-size-factors) for additional details.
- `--normalized_counts_output_prefix` - Flag to specify the location of the normalized counts output files
- `--dge_output_prefix` - Flag to specify the location of the differential gene expression out files
- `--extended_table_output_prefix` - Flag to specify the extended table output files prefix (Note: these files are used internally by GeneLab for the visualization portal)
- `--extended_table_output_suffix` - Flag to specify the extended table output files suffix (Note: these files are used internally by GeneLab for the visualization portal)

**Input Data:**

- path/to/runsheet.csv (runsheet csv file generated in [Step 9](#9-create-runsheet))
- [organisms.csv](https://github.com/J-81/Nextflow_RCP/blob/dev_rc1.0.6/assets/organisms.csv) (csv file containing short name, species name, taxon ID, and annotation db object of model organisms hosted on GeneLab)
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))

**Output Data:**

- RSEM_Unnormalized_Counts.csv\# (table containing raw RSEM gene counts for each sample)
- ERCC_SampleTable.csv\# (table containing samples and their respective groups)
- ERCC_rawCounts_unfiltered.csv (table containing raw ERCC unfiltered counts)
- ERCC_rawCounts_filtered.csv (ERCC counts table after removing ERCC genes with low counts)
- ERCC_Normalized_Counts.csv\# (table containing ERCC-normalized gene counts for each sample)
- visualization_output_table_ERCCnorm.csv (file used to generate GeneLab DGE visualizations for ERCC-normalized data)
- visualization_PCA_table_ERCCnorm.csv (file used to generate GeneLab PCA plots for ERCC-normalized data)
- ERCCnorm_differential_expression.csv\# (table containing ERCC-normalized counts for each sample, group statistics, DESeq2 DGE results for each pairwise comparison, and gene annotations)
- ERCCnorm_contrasts.csv\# (table containing all pairwise comparisons for samples containing ERCC spike-in)

> Note: RNAseq processed data interactive tables and plots are found in the [GLDS visualization portal](https://visualization.genelab.nasa.gov/data/studies).

<br>

---

## 10. Evaluate ERCC Spike-In Data 

<br>

### 10a. Evaluate ERCC Count Data in Python

```python
### Setting up the notebook

# import python packages
import pandas as pd
pd.set_option('mode.chained_assignment', None) # suppress chained indexing warnings
import numpy as np
from urllib.request import urlopen, quote, urlretrieve
from json import loads
from re import search
import zipfile
import seaborn as sns
from scipy.stats import linregress
import matplotlib.pyplot as plt


### Use GeneLab API to locate metadata

GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
GLDS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/data/"
FILELISTINGS_URL_PREFIX = GENELAB_ROOT + "/genelab/data/study/filelistings/"
ISA_ZIP_REGEX = r'.*_metadata_.*[_-]ISA\.zip$'

def read_json(url):
    with urlopen(url) as response:
        return loads(response.read().decode())

def get_isa(accession):
    glds_json = read_json(GLDS_URL_PREFIX + accession)
    try:
        _id = glds_json[0]["_id"]
    except (AssertionError, TypeError, KeyError, IndexError):
        raise ValueError("Malformed JSON?")
    isa_entries = [
        entry for entry in read_json(FILELISTINGS_URL_PREFIX + _id)
        if search(ISA_ZIP_REGEX, entry["file_name"])
    ]
    if len(isa_entries) == 0:
        raise ValueError("Unexpected: no ISAs found")
    elif len(isa_entries) > 1:
        raise ValueError("Unexpected: multiple files match the ISA regex")
    else:
        entry = isa_entries[0]
        version = entry["version"]
        url = GENELAB_ROOT + entry["remote_url"] + "?version={}".format(version)
        alt_url = (
            GENELAB_ROOT + "/genelab/static/media/dataset/" +
            quote(entry["file_name"]) + "?version={}".format(version)
        )
        return entry["file_name"], version, url, alt_url


### Get and parse data and metadata

# Get and unzip ISA.zip to extract metadata.

accession = 'GLDS-NNN' # Replace Ns with GLDS number
isaurl = get_isa(accession)[3]
filehandle, _ = urlretrieve(isaurl)
zip_file_object = zipfile.ZipFile(filehandle, 'r')
zip_file_object.namelist() # Print contents of zip file. Pick relevant one from list

# There are datasets that have multiple assays (including microarray), so the RNAseq ISA files from the above output must be selected. 
# Txt files outputted above are indexed as 0, 1, 2, etc. Fill in the indexed number corresponding to the sample (s_\*txt) and assay files for RNAseq (a_\*_(RNA-Seq).txt) in the code block below.

# Extract metadata from the sample file (s_\*txt)
sample_file = zip_file_object.namelist()[1] # replace [1] with index corresponding to the (s_\*txt) file
file = zip_file_object.open(sample_file)
sample_table = pd.read_csv(zip_file_object.open(sample_file), sep='\t')

# Extract metadata from the assay (a_\*_(RNA-Seq).txt) file
assay_file = zip_file_object.namelist()[0] # replace [0] with index corresponding to the (a_\*_(RNA-Seq).txt) file
file = zip_file_object.open(assay_file)
assay_table = pd.read_csv(zip_file_object.open(assay_file), sep='\t')

# Check the sample table
pd.set_option('max_columns', None)
print(sample_table.head(n=3))

# Check the assay table
pd.set_option('max_columns', None)
assay_table.head(n=3)



# Get raw counts table

raw_counts_table = pd.read_csv('/path/to/RSEM_Unnormalized_Counts.csv', index_col=0) 
raw_counts_table.index.rename('Gene_ID', inplace=True)
print(raw_counts_table.head(n=3))

raw_counts_transcripts = raw_counts_table[raw_counts_table.index.str.contains('^ENSMUSG')] # Change according to organism of interest
raw_counts_transcripts = raw_counts_transcripts.sort_values(by=list(raw_counts_transcripts), ascending=False)
print(raw_counts_transcripts)


# Get ERCC counts

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts.reset_index(inplace=True)
ercc_counts = ercc_counts.rename(columns={'Gene_ID':'ERCC ID'})
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())

# Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
filehandle, _ = urlretrieve(ercc_url)
ercc_table = pd.read_csv(filehandle, '\t')
print(ercc_table.head(n=3))



### Calculate the number of ERCC genes detected in each of the 4 (A, B, C and D) groups for each sample

# Extract ERCC counts and calculate the log(2)

meltERCC = ercc_counts.melt(id_vars=['ERCC ID'])
meltERCC['log2 Count'] = meltERCC['value']+1
meltERCC['log2 Count'] = np.log2(meltERCC['log2 Count'])
meltERCC = meltERCC.rename(columns={'variable':'Sample Name', 'value':'Count'})
print(meltERCC.head(n=3))

# Build Mix dictionary to link sample name to mix added and read depth using the assay table

mix_dict = assay_table.filter(['Sample Name','Parameter Value[Spike-in Mix Number]', 
                       'Parameter Value[Read Depth]'])
mix_dict = mix_dict.rename(columns={'Parameter Value[Spike-in Mix Number]':'Mix',
                                    'Parameter Value[Read Depth]':
                                    'Total Reads'})
print(mix_dict.head(n=3))


# Make combined ercc counts and assay table

merged_ercc = meltERCC.merge(mix_dict, on='Sample Name')
print(merged_ercc)

# Read ERCC info including concentrations from merged_ercc table

groupA = ercc_table.loc[ercc_table['subgroup'] == 'A']['ERCC ID']
groupB = ercc_table.loc[ercc_table['subgroup'] == 'B']['ERCC ID']
groupC = ercc_table.loc[ercc_table['subgroup'] == 'C']['ERCC ID']
groupD = ercc_table.loc[ercc_table['subgroup'] == 'D']['ERCC ID']


# Make a dictionary for ERCC groups

group_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['subgroup']))


# Calculate ERCC counts per million and log(2) counts per million

merged_ercc['Count per million'] = merged_ercc['Count'] / (merged_ercc['Total Reads'] / 1000000.0)
merged_ercc['log2 Count per million'] = np.log2(merged_ercc['Count per million']+1)


# Add ERCC group column

merged_ercc['ERCC group'] = merged_ercc['ERCC ID'].map(group_dict)
merged_ercc = merged_ercc.sort_values(by=['Mix'], ascending=True)
print(merged_ercc)


### Filter and calculate mean counts per million of Mix1 and Mix2 spiked samples in each of the 4 groups

# Filter Mix1 CPM and Mix2 CPM in group A 

Adf = merged_ercc.loc[merged_ercc['ERCC group'] == 'A']
Amix1df = Adf.loc[Adf['Mix']=='Mix 1']
Amix1df['Mix1 CPM'] = Amix1df[Amix1df['Count per million'] > 0]['Count per million'].dropna()
Amix1df = Amix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Amix1df = Amix1df.to_frame()
Amix2df = Adf.loc[Adf['Mix']=='Mix 2']
Amix2df['Mix2 CPM'] = Amix2df[Amix2df['Count per million'] > 0]['Count per million'].dropna()
Amix2df = Amix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Amix2df = Amix2df.to_frame()

adf = Amix1df.merge(Amix2df, on='ERCC ID', suffixes=('', '_2'))
adf = adf.reset_index()
adf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (adf['Avg Mix1 CPM'] / adf['Avg Mix2 CPM'])


# Filter Mix1 CPM and Mix2 CPM in group B

Bdf = merged_ercc.loc[merged_ercc['ERCC group'] == 'B']
Bmix1df = Bdf.loc[Bdf['Mix']=='Mix 1']
Bmix1df['Mix1 CPM'] = Bmix1df[Bmix1df['Count per million'] > 0]['Count per million'].dropna()
Bmix1df = Bmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Bmix1df = Bmix1df.to_frame()
Bmix2df = Bdf.loc[Bdf['Mix']=='Mix 2']
Bmix2df['Mix2 CPM'] = Bmix2df[Bmix2df['Count per million'] > 0]['Count per million'].dropna()
Bmix2df = Bmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Bmix2df = Bmix2df.to_frame()

bdf = Bmix1df.merge(Bmix2df, on='ERCC ID')
bdf = bdf.reset_index()
bdf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (bdf['Avg Mix1 CPM'] / bdf['Avg Mix2 CPM'])


# Filter Mix1 CPM and Mix2 CPM in group C

Cdf = merged_ercc.loc[merged_ercc['ERCC group'] == 'C']
Cmix1df = Cdf.loc[Cdf['Mix']=='Mix 1']
Cmix1df['Mix1 CPM'] = Cmix1df[Cmix1df['Count per million'] > 0]['Count per million'].dropna()
Cmix1df = Cmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Cmix1df = Cmix1df.to_frame()
Cmix2df = Cdf.loc[Cdf['Mix']=='Mix 2']
Cmix2df['Mix2 CPM'] = Cmix2df[Cmix2df['Count per million'] > 0]['Count per million'].dropna()
Cmix2df = Cmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Cmix2df = Cmix2df.to_frame()

cdf = Cmix1df.merge(Cmix2df, on='ERCC ID')
cdf = cdf.reset_index()
cdf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (cdf['Avg Mix1 CPM'] / cdf['Avg Mix2 CPM'])


# Filter Mix1 CPM and Mix2 CPM in group D

Ddf = merged_ercc.loc[merged_ercc['ERCC group'] == 'D']
Dmix1df = Ddf.loc[Ddf['Mix']=='Mix 1']
Dmix1df['Mix1 CPM'] = Dmix1df[Dmix1df['Count per million'] > 0]['Count per million'].dropna()
Dmix1df = Dmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Dmix1df = Dmix1df.to_frame()
Dmix2df = Ddf.loc[Ddf['Mix']=='Mix 2']
Dmix2df['Mix2 CPM'] = Dmix2df[Dmix2df['Count per million'] > 0]['Count per million'].dropna()
Dmix2df = Dmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Dmix2df = Dmix2df.to_frame()

ddf = Dmix1df.merge(Dmix2df, on='ERCC ID')
ddf = ddf.reset_index()
ddf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (ddf['Avg Mix1 CPM'] / ddf['Avg Mix2 CPM'])


##### Multi-sample ERCC analyses

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group A in Mix 1 and Mix 2 spiked samples

a = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupA, hue="Mix",data=merged_ercc[merged_ercc['ERCC ID'].isin(groupA)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
a.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 4")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group A ERCC genes (for group A we expect Mix 1 CPM / Mix 2 CPM = 4)

a1 = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=adf, kind="bar", height=5, aspect=1, linewidth=0.5)
a1.set_xticklabels(rotation=90)
plt.title("ERCC Group A")
a1.set(ylim=(0, 5))
print('Number of ERCC detected in group A (out of 23) =', adf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())


### Create box and whisker plots of the log(2) CPM for each ERCC detected in group B in Mix 1 and Mix 2 spiked samples

b = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupB, hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupB)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
b.set_xticklabels(rotation=90)
plt.text(23,2,"Mix1/ Mix2 = 1")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group B ERCC genes (for group B we expect Mix 1 CPM / Mix 2 CPM = 1)

b = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=bdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
b.set_xticklabels(rotation=90)
plt.title("ERCC Group B")
b.set(ylim=(0, 2))
print('Number of ERCC detected in group B (out of 23) =', bdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())


### Create box and whisker plots of the log(2) CPM for each ERCC detected in group C in Mix 1 and Mix 2 spiked samples

c = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupC, hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupC)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
c.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.67")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group C ERCC genes (for group C we expect Mix 1 CPM / Mix 2 CPM = 0.67)

c = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=cdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
c.set_xticklabels(rotation=90)
plt.title("ERCC Group C")
c.set(ylim=(0, 2.6))
print('Number of ERCC detected in group C (out of 23) =', cdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())


### Create box and whisker plots of the log(2) CPM for each ERCC detected in group D in Mix 1 and Mix 2 spiked samples

d = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupD, hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupD)], col="ERCC group", kind="box", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
d.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.5")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group D ERCC genes (for group D we expect Mix 1 CPM / Mix 2 CPM = 0.5)

d = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=ddf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
d.set_xticklabels(rotation=90)
plt.title("ERCC Group D")
d.set(ylim=(0, 2))
print('Number of ERCC detected in group D (out of 23) =', ddf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())



##### Individual sample ERCC analyses

#Calculate and plot ERCC metrics from individual samples, including limit of detection, dynamic range, and R^2 of counts vs. concentration.
print(ercc_table.head(n=3))

# Make a dictionary for ERCC concentrations for each mix

mix1_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 1 (attomoles/ul)']))
mix2_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 2 (attomoles/ul)']))

# Check assay_table header to identify the 'Sample Name' column and the column title indicating the 'Spike-in Mix Nmber' if it's indicated in the metadata.

pd.set_option('max_columns', None)
print(assay_table.head(n=3))

# Get samples that use mix 1 and mix 2

mix1_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 1']['Sample Name']
mix2_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 2']['Sample Name']

# Get ERCC counts for all samples

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())

# Get ERCC counts for Mix 1 spiked samples

ercc_counts_mix_1 = ercc_counts[mix1_samples]
ercc_counts_mix_1['ERCC conc (attomoles/ul)'] = ercc_counts_mix_1.index.map(mix1_conc_dict)
print(ercc_counts_mix_1.head(n=3))

# Get ERCC counts for Mix 2 spiked samples

ercc_counts_mix_2 = ercc_counts[mix2_samples]
ercc_counts_mix_2['ERCC conc (attomoles/ul)'] = ercc_counts_mix_2.index.map(mix2_conc_dict)
print(ercc_counts_mix_2.head(n=3))


# Create a scatter plot of log(2) ERCC counts versus log(2) ERCC concentration for each sample

columns_mix_1 = ercc_counts_mix_1.columns.drop(['ERCC conc (attomoles/ul)'])
columns_mix_2 = ercc_counts_mix_2.columns.drop(['ERCC conc (attomoles/ul)'])
all_columns = columns_mix_1.to_list() + columns_mix_2.to_list()
total_columns = len(columns_mix_1) + len(columns_mix_2) 
side_size = np.int32(np.ceil(np.sqrt(total_columns)))# calculate grid side size. take sqrt of total plots and round up.
fig, axs = plt.subplots(side_size, side_size, figsize=(15,15), sharex='all', sharey='all'); #change figsize x,y labels if needed.
fig.tight_layout(pad=1, w_pad=2.5, h_pad=3.5)

counter = 0
for ax in axs.flat:
    
    if(counter < len(columns_mix_1)):
      ax.scatter(x=np.log2(ercc_counts_mix_1['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_1[all_columns[counter]]+1), s=7);
      ax.set_title(all_columns[counter][-45:], fontsize=9);
      ax.set_xlabel('log2 ERCC conc (attomoles/ ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=8, labelleft=True, labelbottom=True);
      
    elif(counter >= len(columns_mix_1) and counter < total_columns):
      ax.scatter(x=np.log2(ercc_counts_mix_2['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_2[all_columns[counter]]+1), s=7);
      ax.set_title(all_columns[counter][-45:], fontsize=9);
      ax.set_xlabel('log2 ERCC conc (attomoles/ ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=8, labelleft=True, labelbottom=True);
       
    else:
      pass

    counter = counter + 1


## Calculate and plot linear regression of log(2) ERCC counts versus log(2) ERCC concentration for each sample

# Filter counts > 0

nonzero_counts_list_1 = []
for i in range(0, len(ercc_counts_mix_1.columns)-1):
  counts = ercc_counts_mix_1[columns_mix_1[i]]
  counts.index.rename('Gene_ID', inplace=True)
  countsdf = pd.DataFrame(counts)
  nonzero_counts = countsdf[ercc_counts_mix_1[columns_mix_1[i]] > 0.0]
  nonzero_counts['Conc'] = nonzero_counts.index.map(mix1_conc_dict)
  nonzero_counts.columns = ['Counts','Conc']
  nonzero_counts_sorted = nonzero_counts.sort_values('Conc')
  nonzero_counts_list_1.append(nonzero_counts_sorted)

nonzero_counts_list_2 = []
for i in range(0, len(ercc_counts_mix_2.columns)-1):
  counts = ercc_counts_mix_2[columns_mix_2[i]]
  counts.index.rename('Gene_ID', inplace=True)
  countsdf = pd.DataFrame(counts)
  nonzero_counts = countsdf[ercc_counts_mix_2[columns_mix_2[i]] > 0.0]
  nonzero_counts['Conc'] = nonzero_counts.index.map(mix2_conc_dict)
  nonzero_counts.columns = ['Counts','Conc']
  nonzero_counts_sorted = nonzero_counts.sort_values('Conc')
  nonzero_counts_list_2.append(nonzero_counts_sorted)


# Plot each sample using linear regression of scatter plot with x = log2 Conc and y = log2 Counts.  Return min, max, R^2 and dynamic range (max / min) values.

samples = []
mins = []
maxs = []
dyranges = []
rs = []

fig, axs = plt.subplots(side_size, side_size, figsize=(20,15), sharex='all', sharey='all');
fig.tight_layout(pad=1, w_pad=2.5, h_pad=3.5)

counter = 0
list2counter = 0
for ax in axs.flat:
    
    if(counter < len(columns_mix_1)):

      nonzero_counts = nonzero_counts_list_1[counter]
      xvalues = nonzero_counts['Conc']
      yvalues = nonzero_counts['Counts']

      sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax);
      ax.set_title(all_columns[counter][-47:], fontsize=9);
      ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=8, labelleft=True, labelbottom=True)
      samples.append(all_columns[counter])

      if(len(xvalues) == 0):
        mins.append('NaN')
        maxs.append('NaN')
        dyranges.append('NaN')
        rs.append('NaN')

    
      else:
        min = xvalues[0];
        mins.append(min)
        minimum = f'Min:{min:.1f}';
        max = xvalues[-1];
        maxs.append(max)
        maximum = f'Max:{max:.1f}';
        dynamic_range = max / min;
        dyranges.append(dynamic_range)
        dyn_str = f'Dyn:{dynamic_range:.1f}';

        ax.text(0.02, 0.98, minimum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.88, maximum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.78, dyn_str,verticalalignment='top',
                horizontalalignment='left',transform=ax.transAxes,
                color='black', fontsize=10);
      
        if(len(xvalues) == 1):
          rs.append('NaN')

        else:
          slope, intercept, r, p, se = linregress(np.log2(xvalues), y=np.log2(yvalues))
          r_str = f'R:{r:.2f}'
          rs.append(r)

          ax.text(0.02, 0.68, r_str, verticalalignment='top',
                  horizontalalignment='left',transform=ax.transAxes,
                  color='black', fontsize=10);
    
    elif(counter >= len(columns_mix_1) and counter < total_columns):
      
      nonzero_counts = nonzero_counts_list_2[list2counter]
      xvalues = nonzero_counts['Conc']
      yvalues = nonzero_counts['Counts']

      sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax);
      ax.set_title(all_columns[counter][-47:], fontsize=9);
      ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=8, labelleft=True, labelbottom=True);
      samples.append(all_columns[counter])


      if(len(xvalues) == 0):
        mins.append('NaN')
        maxs.append('NaN')
        dyranges.append('NaN')
        rs.append('NaN')
    
      else:
        min = xvalues[0];
        mins.append(min)
        minimum = f'Min:{min:.1f}';
        max = xvalues[-1];
        maxs.append(max)
        maximum = f'Max:{max:.1f}';
        dynamic_range = max / min;
        dyranges.append(dynamic_range)
        dyn_str = f'Dyn:{dynamic_range:.1f}';

        ax.text(0.02, 0.98, minimum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.88, maximum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.78, dyn_str,verticalalignment='top',
                horizontalalignment='left',transform=ax.transAxes,
                color='black', fontsize=10);
      
        if(len(xvalues) == 1):
          rs.append('NaN')
          
        else:
          slope, intercept, r, p, se = linregress(np.log2(xvalues), y=np.log2(yvalues));
          r_str = f'R:{r:.2f}';
          rs.append(r)

          ax.text(0.02, 0.68, r_str, verticalalignment='top',
                  horizontalalignment='left',transform=ax.transAxes,
                  color='black', fontsize=10);

      list2counter = list2counter + 1
    
    else:
      pass

    counter = counter + 1


# Create directory for saved files

import os
os.makedirs(name="ERCC_analysis", exist_ok=True)

# Print tables containing the dynamic range and R^2 values for each sample.
# Remember to change file names to GLDS# analyzing

stats = pd.DataFrame(list(zip(samples, mins, maxs, dyranges, rs)))
stats.columns = ['Samples', 'Min', 'Max', 'Dynamic range', 'R']
stats.to_csv('ERCC_analysis/ERCC_stats_GLDS-NNN.csv', index = False) 
stats.filter(items = ['Samples', 'Dynamic range']).to_csv('ERCC_analysis/ERCC_dynrange_GLDS-NNN_mqc.csv', index = False) 
stats.filter(items = ['Samples', 'R']).to_csv('ERCC_analysis/ERCC_rsq_GLDS-NNN_mqc.csv', index = False) 



### Generate data and metadata files needed for ERCC DESeq2 analysis

# ERCC Mix 1 and Mix 2 are distributed so that half the samples receive Mix 1 spike-in and half receive Mix 2 spike-in. Transcripts in Mix 1 and Mix 2 are present at a known ratio, so we can determine how well these patterns are revealed in the dataset.

# Get sample table

combined = sample_table.merge(assay_table, on='Sample Name')
combined = combined.set_index(combined['Sample Name'])
pd.set_option('max_columns', None)
print(combined)

# Create metadata table containing samples and their respective ERCC spike-in Mix number
# Sometimes Number in [Spike-in Mix Number] is spelled 'number' and this could cause error in mismatch search 

ERCCmetadata = combined[['Parameter Value[Spike-in Mix Number]']]
ERCCmetadata.index = ERCCmetadata.index.str.replace('-','_')
ERCCmetadata.columns = ['Mix']
#ERCCmetadata = ERCCmetadata.rename(columns={'Parameter Value[Spike-in Mix Number]':'Mix'})
print(ERCCmetadata)

# Export ERCC sample metadata

ERCCmetadata.to_csv('ERCC_analysis/ERCCmetadata.csv') 

# Export ERCC count data

ercc_counts.columns = ercc_counts.columns.str.replace('-','_')
ERCCcounts = ercc_counts.loc[:,ERCCmetadata.index]
ERCCcounts.head()

ERCCcounts.to_csv('ERCC_analysis/ERCCcounts.csv') 
```


**Input Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, used to define sample groups - the \*ISA.zip file is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files')
- RSEM_Unnormalized_Counts.csv (RSEM counts table, output from [Step 9a](#9a-for-datasets-with-ercc-spike-in))

**Output Data:**

- ERCC_analysis/ERCC_stats_GLDS-*.csv (Samplewise counts statistics table containing 'Min', 'Max', 'Dynamic range', 'R')
- ERCC_analysis/ERCC_dynrange_GLDS-*.csv (Samplewise counts statistics subset table containing 'Dynamic range')
- ERCC_analysis/ERCC_rsq_GLDS-*.csv (Samplewise counts statistics subset table containing 'R')
- ERCC_analysis/ERCCmetadata.csv (Samplewise metadata table inlcuding ERCC mix number)
- ERCC_analysis/ERCCcounts.csv (Samplewise ERCC counts table)

<br>

### 10b. Perform DESeq2 Analysis of ERCC Counts in R

```R

## Install R packages if not already installed

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

## Import DESeq2 library

library("DESeq2")

## Import and format ERCC count data and metadata

cts <- as.matrix(read.csv('ERCC_analysis/ERCCcounts.csv',sep=",",row.names="Gene_ID")) #INPUT
coldata <- read.csv('ERCC_analysis/ERCCmetadata.csv', row.names=1) #INPUT

coldata$Mix <- factor(coldata$Mix)
all(rownames(coldata) == colnames(cts))


## Make DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Mix)
dds


## Filter out ERCC genes with counts of less than 10 in all samples #####

keepGenes <- rowSums(counts(dds)) > 10
dds <- dds[keepGenes,]

dds


## Run DESeq2 analysis and calculate results

dds <- DESeq(dds)
res <- results(dds, contrast=c("Mix","Mix 1","Mix 2"))
res

## Export DESeq2 results table and normalized ERCC counts table

write.csv(res, 'ERCC_analysis/ERCC_DESeq2.csv') #OUTPUT
normcounts = counts(dds, normalized=TRUE)
write.csv(normcounts, 'ERCC_analysis/ERCC_normcounts.csv') #OUTPUT
```

**Input Data:**

- ERCC_analysis/ERCCmetadata.csv (Samplewise metadata table inlcuding ERCC mix number, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))
- ERCC_analysis/ERCCcounts.csv (Samplewise ERCC counts table, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))

**Output Data:**

- ERCC_analysis/ERCC_DESeq2.csv (DESeq2 results table)
- ERCC_analysis/ERCC_normcounts.csv (Normalized ERCC Counts table)

<br>

### 10c. Analyze ERCC DESeq2 Results in Python

```python

# Import python packages

import pandas as pd
from urllib.request import urlopen, quote, urlretrieve
import seaborn as sns
import matplotlib.pyplot as plt


# Import ERCC DESeq2 results

deseq2out = pd.read_csv('ERCC_analysis/ERCC_DESeq2.csv', index_col=0) # INPUT
#deseq2out.index = deseq2out.index.str.replace('_','-')
deseq2out.rename(columns ={'baseMean' : 'meanNormCounts'}, inplace = True)
print(deseq2out.head())


# Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
filehandle, _ = urlretrieve(ercc_url)
ercc_table = pd.read_csv(filehandle, '\t', index_col='ERCC ID')
print(ercc_table.head(n=3))


# Combine ERCC DESeq2 results and ercc_table

combined = deseq2out.merge(ercc_table, left_index=True, right_index=True)
print(combined.head())


# Filter p-value and adj. p-value cutoff at 10^-3

combined['cleaned_padj'] = combined['padj']
combined.loc[(combined.cleaned_padj < 0.001),'cleaned_padj']=0.001

combined['cleaned_pvalue'] = combined['pvalue']
combined.loc[(combined.cleaned_pvalue < 0.001),'cleaned_pvalue']=0.001

print(combined.head())


# Export the filtered combined ERCC DESeq2 results and ercc_table
# Remember to change file name to GLDS# analyzing

combined.filter(items = ['ERCC ID', 'meanNormCounts', 'cleaned_pvalue','cleaned_padj']).to_csv('ERCC_analysis/ERCC_lodr_GLDS-NNN_mqc.csv') 


# Plot p-value vs. mean normalized ERCC counts

fig, ax = plt.subplots(figsize=(10, 7))

sns.scatterplot(data=combined, x="meanNormCounts", y="cleaned_pvalue",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

sns.lineplot(data=combined, x="meanNormCounts", y="cleaned_pvalue",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

#g.set_xscale("log", base=2)
ax.set_xscale("linear");
ax.set_yscale("log");


# Plot Adjp-value vs. mean normalized ERCC counts

fig, ax = plt.subplots(figsize=(10, 7))

sns.scatterplot(data=combined, x="meanNormCounts", y="cleaned_padj",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

sns.lineplot(data=combined, x="meanNormCounts", y="cleaned_padj",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

#g.set_xscale("log", base=2)
ax.set_xscale("linear");
ax.set_yscale("log");
```

**Input Data:**

- ERCC_analysis/ERCC_DESeq2.csv (ERCC DESeq2 results table, output from [Step 10b](#10b-perform-deseq2-analysis-of-ercc-counts-in-r))

**Output Data:**

- ERCC_analysis/ERCC_lodr_*.csv (ERCC Gene Table including mean counts, adjusted p-value and p-value, and filtered to genes with both adj. p-value and p-value < 0.001)

> All steps of the ERCC Spike-In Data Analysis are performed in a Jupyter Notebook (JN) and the completed JN is exported as an html file and published in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) for the respective dataset.
