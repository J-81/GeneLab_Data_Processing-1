# NF_RCP-F Workflow Information and Usage Instructions

## General Workflow Info

### Implementation Tools

The current GeneLab RNAseq consensus processing pipeline (RCP), [GL-DPPD-7101-F](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in containers. This workflow is run using the CLI of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.

### Workflow & Subworkflows

---

- **Click image to expand**

<p align="center">
<a href="../../images/NF_RCP-F_rnaseq_workflow.png"><img src="../../images/NF_RCP-F_rnaseq_workflow.png"></a>
</p>

---
The NF_RCP-F workflow is composed of three subworkflows as shown in the image above.
Below is a description of each subworkflow and the additional output files generated that are not already indicated in the [GL-DPPD-7101-F pipeline 
document](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md):

1. **Analysis Staging Subworkflow**

   - Description:
     - This subworkflow extracts the processing parameters (e.g. organism, library layout) from the GLDS ISA archive as well as retrieves the raw reads files hosted on the GeneLab Data Repository.

2. **RNASeq Consensus Pipeline Subworkflow**

   - Description:
     - This subworkflow uses the staged raw data and processing parameters to generate processed data.

3. **V&V Pipeline Subworkflow**

   - Description:
     - This subworkflow performs validation and verification on the raw and processed files.  It performs a series of checks and flags the results to a series of log files. The following flag levels are found in these logs:

---
|Flag Codes|Flag Name|Interpretation|
|:---------|:--------|:-------------|
| 20-29 | GREEN | Indicates the check passed all validation conditions |
| 30-39 | YELLOW | Indicates the check was flagged for minor issues (e.g. slight outliers) |
| 50-59 | RED | Indicates the check was flagged for moderate issues (e.g. major outliers) |
| 80-89 | HALT | Indicates the check was flagged for severe issues that trigger a processing halt (e.g. missing data) |

## Utilizing the Workflow


### 1. Install Nextflow and Singularity

#### 1a. Installing Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the Nextflow documentation [page](https://www.nextflow.io/docs/latest/getstarted.html).

#### 1b. Installing Singularity

Singularity is a container platform that allows usage of containerized software. This enables the workflow to retrieve and use all software required during the processing pipeline without needing to directly install the processing software directly to the user system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

### 2. Download the Workflow Files

All files required for utilizing the NF_RCP-E GeneLab workflow for processing RNASeq data are in the [workflow_code](workflow_code) directory. To get a 
copy of latest NF_RCP-E version on to your system, copy the github web address of the [latest NF_RCP-E version](workflow_code/NF_RCP-E_1.0.0) 
sub-directory under the workflow_code 
directory, then paste it into [GitZip here](http://kinolien.github.io/gitzip/), and click download:

<p align="center">
<a href="../../images/NF_RCP-F_gitzip_rnaseq.png"><img src="../../images/NF_RCP-F_gitzip_rnaseq.png"></a>
</p>

### 3. Setup Execution Permission for Bin Scripts

Once you've downloaded the workflow template, you need to set the execution permission for the scripts in the bin folder.  The scripts may be made executable using the following command inside the unzipped workflow_code directory.

```bash
chmod -R u+x bin
```

### 4. Run the Workflow

#### Approach 1: Running the workflow with automatic retrieval of Ensembl reference fasta and gtf

Here is one example command of how to run the workflow in using Approach 1.  Note: main.nf is a file located in the workflow_code directory.

> **Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --ensemblVersion) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument**  

``` text
Usage example 1:
   Fetches ensembl reference files via ftp and GeneLab raw data via https before running processing pipeline
   > nextflow run ./main.nf --gldsAccession GLDS-194 --ensemblVersion 96

Usage example 2:
   Fetches GeneLab raw data via https before running processing pipeline using supplied local reference fasta and gtf files.
   Note: ensemblVersion and ref_source are used here to label subdirectories for derived reference files.
   > nextflow run ./main.nf --gldsAccession GLDS-194 --ensemblVersion 96 --ref_source <reference_label>  --ref_fasta </path/to/fasta> --ref_gtf </path/to/gtf>

required arguments:
  --gldsAccession GLDS-000
                        the GLDS accession id to process through the RNASeq consensus Pipeline.
  --ensemblVersion n    the ensembl Version to use for the reference genome.
optional arguments:
  --help                show this help message and exit
  --skipVV              skip automated V&V processes. Default: false
  --outputDir           directory to save staged raw files and processed files. Default: <launch directory>
  --limitSamplesTo n    limit the number of samples staged to a number.
  --genomeSubsample n   subsamples genome fasta and gtf files to the supplied chromosome.
  --truncateTo n        limit number of reads downloaded and processed to *n* reads , for paired end limits number of reverse and forward read files to *n* reads each.
  --force_single_end    forces analysis to use single end processing.  For paired end datasets, this means only R1 is used.  For single end studies, this should have no effect.
  --stageLocal          download the raw reads files for the supplied GLDS accession id.  Set to false to disable raw read download and processing.  Default: true
  --ref_order           specifies the reference to use from ensembl.  Allowed values:  ['toplevel','primary_assemblyELSEtoplevel']. 'toplevel' : use toplevel.  'primary_assemblyELSEtoplevel' : use primary assembly, but use toplevel if primary assembly doesn't exist. Default: 'primary_assemblyELSEtoplevel'
  --ref_fasta           specifies a reference fasta from a local path. This an is an alternative approach from the automatic retrieval of reference files from ensembl
  --ref_gtf             specifies a reference gtf from a local path. This an is an alternative approach from the automatic retrieval of reference files from ensembl
  --referenceStorePath  specifies the directory where fetched reference files are downloaded to
  --derivedStorePath    specifies the directory where derivative reference files are saved. Examples of such files in this pipeline included BED and PRED files generated from the reference gtf
  --ref_source          a string to label subdirectories in 'StorePath' paths. Examples include 'ensembl' or 'ensembl_plants'.
  -stub-run             runs the workflow forcing 'unstranded' RSEM settings and using dummy gene counts in the differential gene expression (DGE) analysis. Useful when combined with the --truncateTo parameter this often leads to low gene counts and errors in the DGE analysis
```

See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details common to all nextflow workflows.

---

### 5. Additional Output Files

The output from the Analysis Staging subworkflow, V&V Pipeline subworkflow, and Nextflow specific logs are described here.
> Note: The outputs from version F of the RNASeq Consensus Pipeline are documented in the current processing protocol, 
[GL-DPPD-7101-F.md](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md).

1. Analysis Staging Subworkflow

   - Output:
     - \*_bulkRNASeq_v1_runsheet.csv (a table that include processing parameters and raw reads files location)
     - \*-ISA.zip (the ISA archive fetched from the GeneLab Data Repository)
     - \*_metadata_table.txt (a table that includes additional information about the GLDS entry, not used for processing)

1. V&V Pipeline Subworkflow

   - Output:
     - VV_Logs/VV_log_final.tsv (A tab-separated values file that includes all V&V flags levels logged)
     - VV_Logs/VV_log_final_only_issues.tsv (A tab-separated values file that includes V&V flags levels logged with maximum flag codes greater than 20)
     - VV_Logs/VV_log_verbose_through_VV_RAW_READS.tsv (A tab-separated values file that includes all V&V flags levels logged, generated after RAW_READS)
     - VV_Logs/VV_log_verbose_through_VV_TRIMMED_READS.tsv (A tab-separated values file that includes all V&V flags levels logged, generated after TRIMMED_READS)
     - VV_Logs/VV_log_verbose_through_VV_STAR_ALIGNMENTS.tsv (A tab-separated values file that includes all V&V flags levels logged, generated after STAR_ALIGNMENTS)
     - VV_Logs/VV_log_verbose_through_VV_RSEQC.tsv (A tab-separated values file that includes all V&V flags levels logged, generated after RSEQC)
     - VV_Logs/VV_log_verbose_through_VV_RSEM_COUNTS.tsv (A tab-separated values file that includes all V&V flags levels logged, generated after RSEM_COUNTS)

---

### 6. Known Issues to Look Out For

#### Truncated Raw Read Files

- This is a known issue for Nextflow file staging from URL.  
- If the Nextflow process is forcefully interrupted while staging a file (most notably raw read files for this pipeline), the truncated file will **NOT** be re-downloaded, resulting in a pipeline trying to process with the truncated file.
- This is most commonly manifests as an unexpected error related to truncation occuring for processes that use the raw read files.
- The advised workaround is to purge the staged files, located in your Nextflow "work" directory under the "stage" sub-directory, and relaunch the pipeline.
