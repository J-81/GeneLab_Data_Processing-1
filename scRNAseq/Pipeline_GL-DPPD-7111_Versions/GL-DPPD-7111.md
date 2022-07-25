# GeneLab bioinformatics processing pipeline for single cell RNA-sequencing data

> This page holds an overview and instructions for how GeneLab processes single cell RNA-sequencing (scRNAseq) datasets. Exact processing commands and GL-DPPD-7111 version used for specific datasets are available in the 
[GLDS_Processing_Scripts](../GLDS_Processing_Scripts) directory and processed data output files are provided in the [GeneLab Data Systems 
(GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** July XX, 2022  
**Revision:** -  
**Document Number:** GL-DPPD-7111  

**Submitted by:**  
Lauren Sanders (GeneLab Data Processing Team)

**Approved by:**  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Jonathan Galazka (GeneLab Project Scientist)

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [**2. Trim/Filter Raw Data and Trimmed Data QC**](#2-trimfilter-raw-data-and-trimmed-data-qc)
    - [2a. Trim/Filter Raw Data](#2a-trimfilter-raw-data)
    - [2b. Trimmed Data QC](#2b-trimmed-data-qc)
    - [2c. Compile Trimmed Data QC](#2c-compile-trimmed-data-qc)
  - [**3. Build STAR Reference**](#3-build-star-reference)
  - [**4. Align Reads to Reference Genome**](#4-align-reads-to-reference-genome)
    - [4a. Align Reads to Reference Genome with STARsolo](#4a-align-reads-to-reference-genome-with-starsolo)
    - [4b. Compile Alignment Logs](#4b-compile-alignment-logs)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.11.9|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.12|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|3.7|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|0.6.7|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|2.7.10a|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are provided in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) directory.
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
- `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are 
immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

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
- `--genomeSAindexNbases` - length (in bases) of the SA pre-indexing string, usually between 10 and 15. Longer strings require more memory but allow for faster searches. This value should be scaled down for smaller genomes (like bacteria) to min(14, 
log2(GenomeLength)/2 - 1). For example, for a 1 megaBase genome this value would be 9.
- `--genomeDir` - specifies the path to the directory where the STAR reference will be stored. At least 100GB of available disk space is required for mammalian genomes.
- `--genomeFastaFiles` - specifies one or more fasta file(s) containing the genome reference sequences
- `--sjdbGTFfile` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `--sjdbOverhang` - indicates the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. The length should be one less than the maximum length of the reads.

**Input Data:**

- *.fasta (genome sequence, this scRCP version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)
- *.gtf (genome annotation, this scRCP version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

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

## 4. Align Reads to Reference Genome

<br>

### 4a. Align Reads to Reference Genome with STARsolo

```bash
STAR --runThreadN <NumberOfThreads> \
  --genomeDir /path/to/STAR/genome/directory \
  --soloType CB_UMI_Simple \ # Used for 10X Chromium data 
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloUMIlen 12 \
  --soloCellFilter EmptyDrops_CR <ExpectedCells> 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
  --soloMultiMappers EM \
  --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM \
  --outSAMtype BAM Unsorted \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --readFilesCommand zcat \
  --soloCBwhitelist CellBarcodeWhitelist \
  --outFileNamePrefix /path/to/STAR/output/directory/<sample_id> \
  --readFilesIn /path/to/reads_from_sample \
  /path/to/barcode_reads

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
- `--quantMode` - specifies the type(s) of quantification desired; the `TranscriptomeSAM` option instructs STAR to output a separate sam/bam file containing alignments to the transcriptome and the `GeneCounts` option instructs STAR to output a tab 
delimited file containing the number of reads per gene
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
- *ReadsPerGene.out.tab (tab deliminated file containing STAR read counts per gene with 4 columns that correspond to different strandedness options: column 1 = gene ID, column 2 = counts for unstranded RNAseq, column 3 = counts for 1st read strand 
aligned with RNA, column 4 = counts for 2nd read strand aligned with RNA)
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

   
