# NF_RCP-F Workflow Information and Usage Instructions

## General Workflow Info

### Implementation Tools

The current GeneLab RNAseq consensus processing pipeline (RCP), [GL-DPPD-7101-F](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in containers. This workflow (NF_RCP-F) is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.

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
     - This subworkflow extracts the metadata parameters (e.g. organism, library layout) needed for processing from the GLDS ISA archive and retrieves the raw reads files hosted on the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects).
       > *GLDS ISA archive*: ISA directory containing Investigation, Study, and Assay (ISA) metadata files for a respective GLDS dataset - the *ISA.zip file is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files' for any GLDS dataset in the GeneLab Data Repository.

2. **RNASeq Consensus Pipeline Subworkflow**

   - Description:
     - This subworkflow uses the staged raw data and metadata parameters from the Analysis Staging Subworkflow to generate processed data using [version F of the GeneLab RCP](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md).

3. **V&V Pipeline Subworkflow**

   - Description:
     - This subworkflow performs validation and verification (V&V) on the raw and processed data files in real-time.  It performs a series of checks on the output files generated and flags the results, using the flag codes indicated in the table below, which outputted as a series of log files. 

       **V&V Flags**:

      | Flag Name             |Interpretation                  |
      |-----------------------|-------------------------|
      | GREEN                 |Indicates the check passed all validation conditions                         |  
      | YELLOW                |Indicates the check was flagged for minor issues (e.g. slight outliers)                         |
      | RED                   |Indicates the check was flagged for moderate issues (e.g. major outliers)                               |
      | HALT                  |Indicates the check was flagged for severe issues that trigger a processing halt (e.g. missing data)                         |

<br>


---

## Utilizing the Workflow


### 1. Install Nextflow and Singularity

#### 1a. Installing Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the Nextflow documentation [page](https://www.nextflow.io/docs/latest/getstarted.html).

#### 1b. Installing Singularity

Singularity is a container platform that allows usage of containerized software. This enables the workflow to retrieve and use all software required during the processing pipeline without needing to directly install the processing software directly to the user system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

### 2. Download the Workflow Files

All files required for utilizing the NF_RCP-F GeneLab workflow for processing RNASeq data are in the [workflow_code](workflow_code) directory. To get a 
copy of latest NF_RCP-F version on to your system, copy the github web address of the [latest NF_RCP-F version](workflow_code/NF_RCP-F_1.0.0), then paste it into [GitZip here](http://kinolien.github.io/gitzip/), and click download:

TODO: Update image when we have the official NASA GitHub link - alternatively create script that can be run to do this automatically
<p align="center">
<a href="../../images/NF_RCP-F_gitzip_rnaseq.png"><img src="../../images/NF_RCP-F_gitzip_rnaseq.png"></a>
</p>

### 3. Setup Execution Permission for Workflow Scripts

Once you've downloaded the NF_RCP-F workflow directory as a zip file, unzip the workflow then `cd` into the NF_RCP-F directory on the CLI. Next, run the following command to set the execution permissions for all scripts in the bin folder:

```bash
chmod -R u+x bin
```

### 4. Run the Workflow

While in the NF_RCP-F workflow directory, you are now able to run the workflow. Below are three examples of how to run the NF_RCP-F workflow:
> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --ensemblVersion) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

**Approach 1: Run the workflow with automatic retrieval of Ensembl reference fasta and gtf files**

```bash
nextflow run ./main.nf --gldsAccession GLDS-194
```

**Approach 2: Run the workflow using local Ensembl reference fasta and gtf files**

```bash
nextflow run ./main.nf --gldsAccession GLDS-194 --ensemblVersion 107 --ref_fasta </path/to/fasta> --ref_gtf </path/to/gtf>
```

**Approach 3: Run the workflow with user-created runsheet**

Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).
```bash
nextflow run ./main.nf --runsheetPath </path/to/runsheet>
```

**Required Arguments:**

* `--gldsAccession GLDS-###` – specifies the GLDS dataset to process through the RCP workflow (replace ### with the GLDS number)
* `--ensemblVersion` - specifies the Ensembl version to use for the reference genome (TODO: There should be a default ensemblVersion that is consistent with the ensembl version used for the RCP version the workflow is running, so this can become an optional argument)
  
  
**Optional Arguments:**

* `--help` – show the NF_RCP-F workflow help menu
* `--skipVV` - skip the automated V&V processes (Default: the automated V&V processes are active) 
* `--outputDir` - specifies the directory to save the raw and processed data files (Default: files are saved in the launch directory)
* `--limitSamplesTo` - specifies the number of samples to process (Default: all samples in the GLDS dataset indicated are processed)
* `--genomeSubsample` - specifies an individual chromosome to use for processing (Default: all chromosomes in the fasta/gtf files are used)
* `--truncateTo` - specifies the number of reads to download and process for each sample (Default: all reads are used)
* `--force_single_end` - forces the analysis to use single end processing; for paired end datasets, this means only R1 is used; for single end datasets, this should have no effect
* `--stageLocal TRUE|FALSE` - TRUE = download the raw reads files for the GLDS dataset indicated, FALSE = disable raw reads download and processing (Default: TRUE)
* `--use_dummy_gene_counts` - Use random gene counts instead of actual RSEM gene counts during DGE. This is useful when using a low number of reads / sumbsampled genomes to prevent DGE errors due to low/zero gene counts.
* `--ref_order toplevel|primary_assemblyELSEtoplevel` - specifies which Ensembl fasta file to use, toplevel = use the toplevel fasta, primary_assemblyELSEtoplevel = use the primary_assembly fasta if available but if not, use the toplevel fasta (Default: primary_assemblyELSEtoplevel)
* `--ref_fasta` - specifices the path to a local fasta file (Default: fasta file is downloaded from Ensembl)
* `--ref_gtf` - specifices the path to a local gtf file (Default: gtf file is downloaded from Ensembl)
* `--referenceStorePath` - specifies the directory to store the Ensembl fasta and gtf files (Default: within the directory structure created by default in the launch directory)
* `--derivedStorePath` - specifies the directory to store the tool-specific indices created during processing (Default: within the directory structure created by default in the launch directory)
* `--ref_source` - specifies the ensembl database source (e.g. ensembl, ensembl_plants, ensembl_bacteria) when local reference files are used, required when the `--ref_fasta` and/or `--ref_gtf` options are specified
* `-stub-run` - forces the workflow to use the RSEM `--strandedness none` setting, used in combination with the `--truncateTo` option to avoid V&V errors when evaluating counts and DGE output files during testing

See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details common to all nextflow workflows.

---

### 5. Additional Output Files

The outputs from the Analysis Staging and V&V Pipeline Subworkflows are described below:
> Note: The outputs from the RNASeq Consensus Pipeline Subworkflow are documented in the [GL-DPPD-7101-F.md](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md) processing protocol.

**Analysis Staging Subworkflow**

   - Output:
     - \*_bulkRNASeq_v1_runsheet.csv (table containing metadata required for processing, including the raw reads files location)
     - \*-ISA.zip (the ISA archive of the GLDS datasets to be processed, downloaded from the GeneLab Data Repository)
     - \*_metadata_table.txt (table that includes additional information about the GLDS dataset, not used for processing)
   
   
**V&V Pipeline Subworkflow**

   - Output:
     - VV_Logs/VV_log_final.tsv (table containing V&V flags for all checks performed)
     - VV_Logs/VV_log_final_only_issues.tsv (table containing V&V flags ONLY for checks that produced a flag level >= 30)
     - VV_Logs/VV_log_verbose_through_VV_RAW_READS.tsv (table containing V&V flags ONLY for raw reads checks)
     - VV_Logs/VV_log_verbose_through_VV_TRIMMED_READS.tsv (table containing V&V flags through trimmed reads checks ONLY)
     - VV_Logs/VV_log_verbose_through_VV_STAR_ALIGNMENTS.tsv (table containing V&V flags through alignment file checks ONLY)
     - VV_Logs/VV_log_verbose_through_VV_RSEQC.tsv (table containing V&V flags through RSeQC file checks ONLY)
     - VV_Logs/VV_log_verbose_through_VV_RSEM_COUNTS.tsv (table containing V&V flags through RSEM raw count file checks ONLY)

Standard Nextflow resource usage logs are also produced as follows:
Further details about these logs can also found within this Nextflow documentation [page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

**Nextflow Resource Usage Logs**

   - Output:
     - Resource_Usage/execution_report_{timestamp}.html (an html report which includes many useful metrics about a workflow execution)
     - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in your pipeline)
     - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

---

### 6. Known Issues to Look Out For

#### Recent Nextflow versions introduced a bug that causes the workflow to fail during retrieval of raw reads files from GeneLab Repository

- [Github Issue Link](https://github.com/nextflow-io/nextflow/issues/2918)
- This will be fixed in an upcoming release of Nextflow. In the meantime, the workflow should work with Nextflow Version 21.10.6, which predates the introduction of the bug.
  - We recommend setting the environment variable 'NXF_VER=21.10.6' to allow Nextflow to automatically update/downgrade to that version on launch.
