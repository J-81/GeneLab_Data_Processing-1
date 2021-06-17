# Workflow information and usage instructions


## General workflow info
The current processing protocol is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments. This workflow is run using the CLI of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow. An introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).  

## Utilizing the workflow

1. [Install conda and Nextflow](#1-install-conda-and-nextflow)  
2. [Download the workflow files](#2-download-the-workflow-files)  
3. [Setup execution permission for bin scripts](#3-setup-execution-permission-for-bin-scripts)  
4. [Run the workflow](#4-run-the-workflow)

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
