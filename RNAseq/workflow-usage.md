# Workflow information and usage instructions


## General workflow info
### Implemenation Tools
The current processing protocol is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments. This workflow is run using the CLI of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow. An introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).  

### Workflow & Subworkflows
---

<p align="center">
<a href="images/rnaseq_pipeline.png"><img src="images/rnaseq_pipeline.png"></a>
</p>

---
This workflow is composed of three subworkflows as shown in the image above.
Below is a description of each subworkflow as well as the output files if not listed in the processing protocol:
1. Analysis Staging Subworkflow
  - Description: 
    - This subworkflow extracts the processing parameters (e.g. organism, library layout) from the GLDS ISA archive as well as retrieves the raw reads files hosted on the GeneLab Data Repository.

2. RNASeq Concensus Pipeline Subworkflow
  - Description:
     - This subworkflow uses the staged raw data and processing parameters to generate processed data.
3. V&V Pipeline Subworkflow 
  - Description:
    - This subworkflow performs validation and verification on the raw and processed files.  It performs a series of checks and flags the results to a series of log files. The following flag levels are found in these logs:
---
| Flag ID | Severity              |
|---------|-----------------------|
| 20      | Info-Only             |
| 30      | Passed-Green          |
| 50      | Warning-Yellow        |
| 60      | Warning-Red           |
| 90      | Issue-Halt_Processing |



## Utilizing the workflow

1. [Install conda and Nextflow](#1-install-conda-and-nextflow)  
2. [Download the workflow files](#2-download-the-workflow-files)  
3. [Setup execution permission for bin scripts](#3-setup-execution-permission-for-bin-scripts)  
4. [Run the workflow](#4-run-the-workflow)
5. [Additional Output Files](#5-additional-output-files)



### 1. Install conda and Nextflow
We recommend installing a Miniconda, Python3 version appropriate for your system, as exemplified in [the above link](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  

Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:

```bash
conda install -c bioconda nextflow
nextflow self-update
```

### 2. Download the workflow files
All files required for utilizing the GeneLab workflow for processing RNASeq data are in the [workflow_code](workflow_code) directory. To get a copy of that directory on to your system, copy the github web address of that directory, paste it into [GitZip here](http://kinolien.github.io/gitzip/), and then click download:

<p align="center">
<a href="images/gitzip_rnaseq.png"><img src="images/gitzip_rnaseq.png"></a>
</p>

### 3. Setup execution permission for bin scripts
Once you've downloaded the workflow template, you need to set the execution permission for the scripts in the bin folder.  The scripts may be made executable using the following command inside the unzipped workflow_code directory.

```bash
chmod -R u+x bin
```

### 4. Run the workflow

Here is one example command of how to run the workflow.  Note: main.nf is a file located in the workflow_code directory.

> **Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --ensemblVersion) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument**  


```bash
nextflow run <path/to/main.nf> --gldsAccession=GLDS-194 --ensemblVersion=96 [--outputDir] [--skipVV] [--limitSamplesTo=<n>] [--genomeSubsample=<n>] [--truncateTo=<n>] [--stageLocal]
```

**Required Parameters**
* `--gldsAccession` – specifies the GLDS accession id to process through the RNASeq Concensus Pipeline
* `--ensemblVersion` – indicates the Ensembl Version to use for the reference genome

**Optional Parameters**
* `--outputDir` – specifies the directory to save staged raw files and processed files. Default: '.'

**Debug Parameters**
> Note: These parameters are offered to allow test workflows to run with minimal compute resources or skip portions of the workflow. They should **NEVER** be used to generate processed data intended for real analysis.

* `--skipVV` – skip automated V&V processes. Default: false
* `--limitSamplesTo=<n>` – limit the number of samples staged to a *n* samples
* `--genomeSubsample=<n>` – subsamples genome fasta and gtf files to the supplied *n* chromosome
* `--truncateTo=<n>` – limit number of reads downloaded and processed to *n* reads , for paired end limits number of reverse and forward read files to *n* reads each
* `--stageLocal` – download the raw reads files for the supplied GLDS accession id.  Set to false to disable raw read download and processing.  Default: true
* `-stub-run` – runs the workflow using dummy gene counts in the differential gene expression (DGE) analysis. Useful when combined with the --truncateTo parameter this often leads to low gene counts and errors in the DGE analysis


See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details common to all nextflow workflows.

---

### 5. Additional Output Files

The output from the Analysis Staging subworkflow and V&V Pipeline subworkflow are described here.
Note: the output from the RNASeq Concensus Pipeline are documented in the current processing protocol, (GL-DPPD-7101-C.md)](GL-DPPD-7101-C.md),

1. Analysis Staging Subworkflow
  - Output:
    - \*_runsheet.csv (a table that include processing parameters and raw reads files location)
    - \*-ISA.zip (the ISA archive fetched from the GeneLab Data Repository)
    - \*_metadata_table.txt (a table that includes additional information about the GLDS entry, not used for processing)

2. V&V Pipeline Subworkflow 
  - Output:
    - VV_Log/VV_FULL_OUT.tsv (A tab-separated values file that includes all V&V flags levels logged)
    - VV_Log/only-issues__VV_FULL_OUT.tsv (A tab-separated values file that includes V&V flags levels logged with severities greater than 30)
    - VV_Log/Summary.tsv (A tab-separated values file that summarizes the percent of samples that have Warnings for each step)
    - VV_Log/all-by-sample.txt (A text file that lists, by sample, all flags with severities greater than 30)
    - VV_Log/bySample/{sample_name}__VV_FULL_OUT.tsv (A series of tab-separated values files that subset the full flag log by sample)
    - VV_Log/byStep/{step_name}__VV_FULL_OUT.tsv (A series of tab-separated values files that subset the full flag log by processing step)

