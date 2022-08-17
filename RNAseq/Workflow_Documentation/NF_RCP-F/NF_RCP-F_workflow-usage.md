# NF_RCP-F Workflow Information and Usage Instructions <!-- omit in toc -->

## General Workflow Info <!-- omit in toc -->

### Implementation Tools <!-- omit in toc -->

The current GeneLab RNAseq consensus processing pipeline (RCP), [GL-DPPD-7101-F](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in containers. This workflow (NF_RCP-F) is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.

### Workflow & Subworkflows <!-- omit in toc -->

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
     - This subworkflow performs validation and verification (V&V) on the raw and processed data files in real-time.  It performs a series of checks on the output files generated and flags the results, using the flag codes indicated in the table below, which are outputted as a series of log files.

        **V&V Flags**:

      | Flag Name             |Interpretation                  |
      |-----------------------|-------------------------|
      | GREEN                 |Indicates the check passed all validation conditions                         |  
      | YELLOW                |Indicates the check was flagged for minor issues (e.g. slight outliers)                         |
      | RED                   |Indicates the check was flagged for moderate issues (e.g. major outliers)                               |
      | HALT                  |Indicates the check was flagged for severe issues that trigger a processing halt (e.g. missing data)                         |

<br>


---

## Utilizing the Workflow <!-- omit in toc -->

- [1. Install Nextflow and Singularity](#1-install-nextflow-and-singularity)
  - [1a. Installing Nextflow](#1a-installing-nextflow)
  - [1b. Installing Singularity](#1b-installing-singularity)
- [2. Download the Workflow Files](#2-download-the-workflow-files)
- [3. Fetch Singularity Images](#3-fetch-singularity-images)
- [4. Run the Workflow](#4-run-the-workflow)
  - [**Approach 1**: Run the workflow with automatic retrieval of the same Ensembl reference fasta and gtf files used by GeneLab](#approach-1-run-the-workflow-with-automatic-retrieval-of-the-same-ensembl-reference-fasta-and-gtf-files-used-by-genelab)
  - [**Approach 2**: Run the workflow on GLDS datasets using local Ensembl reference fasta and gtf files](#approach-2-run-the-workflow-on-glds-datasets-using-local-ensembl-reference-fasta-and-gtf-files)
  - [**Approach 3**: Run the workflow with user-created runsheet](#approach-3-run-the-workflow-with-user-created-runsheet)
- [5. Additional Output Files](#5-additional-output-files)
- [6. Known Issues to Look Out For](#6-known-issues-to-look-out-for)
  - [Recent Nextflow versions introduced a bug that causes the workflow to fail during retrieval of files from GeneLab Repository](#recent-nextflow-versions-introduced-a-bug-that-causes-the-workflow-to-fail-during-retrieval-of-files-from-genelab-repository)

### 1. Install Nextflow and Singularity

#### 1a. Installing Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the Nextflow documentation [page](https://www.nextflow.io/docs/latest/getstarted.html).

#### 1b. Installing Singularity

Singularity is a container platform that allows usage of containerized software. This enables the workflow to retrieve and use all software required during the processing pipeline without needing to directly install the processing software directly to the user system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

### 2. Download the Workflow Files
All files required for utilizing the NF_RCP-F GeneLab workflow for processing RNASeq data are in the [workflow_code](workflow_code) directory.

The code is most easily downloaded from the release page.

An example of downloading and unzipping the code is provided below:

```bash
wget https://github.com/asaravia-butler/GeneLab_Data_Processing/releases/download/NF_RCP-F_1.0.0/NF_RCP-F_1.0.0.zip

unzip NF_RCP-F_1.0.0.zip
```

### 3. Fetch Singularity Images

While Nextflow is able to manage fetching Singularity images from url. This feature currently causes frequent intermittent issues ([nextflow-issue #1210](https://github.com/nextflow-io/nextflow/issues/1210)).

To address this issue, Singularity images should be fetched using the following bash script before running the workflow as follows:

```bash
bash NF_RCP-F_1.0.0/bin/prepull_singularity.sh NF_RCP-F_1.0.0/config/software/by_docker_image.config
```

Fetching takes around 20 minutes depending on your network speed.  Once fetched, the images can be reused across workflow runs.

This will create a 'singularity' folder containing the Singularity images.  This folder should be exported as an Nextflow configuration env variable as follows to ensure Nextflow can locate the fetched images.

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

### 4. Run the Workflow

Below are three examples of how to run the NF_RCP-F workflow:
> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --ensemblVersion) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

#### **Approach 1**: Run the workflow with automatic retrieval of the same Ensembl reference fasta and gtf files used by GeneLab

```bash
nextflow run NF_RCP-F_1.0.0/main.nf \
  -profile singularity \
  --gldsAccession GLDS-194 \
  --ensemblVersion 107 \
  --ref_source ensembl
```

**Parameter Definitions:**

- `NF_RCP-F_1.0.0/main.nf` - The nextflow workflow file
- `--gldsAccession` - The GLDS accession ID to process
- `-profile` – Designates the configuration profile(s) to load
- `singularity` - A configuration profile to setup and use singularity for all software
- `--ensemblVersion` – Indicates the Ensembl Release Version used. This is used to organize the storage of reference genome files.
- `107` - The Ensembl Release Version GeneLab uses for Mus Musculus.
  - Other organisms currently use the following versions:
- `--ref_source` – Indicates the reference genome source. This is used to organize the storage of reference genome and indices files.
- `ensembl` - The Ensembl Release Version GeneLab uses for Mus Musculus.
  - Organisms currently supported through the entire workflow use the following versions:

    | Organism                                         | ensemblVersion | ref_source       |
    | :----------------------------------------------- | :------------- | :--------------- |
    | Arabidopsis thaliana                             | 54             | ensembl_plants   |
    | Bacillus subtilis (strain168)                    | 54             | ensembl_bacteria |
    | Caenorhabditis elegans                           | 107            | ensembl          |
    | Danio rerio                                      | 107            | ensembl          |
    | Drosophila melanogaster                          | 107            | ensembl          |
    | Homo sapiens                                     | 107            | ensembl          |
    | Mus musculus                                     | 107            | ensembl          |
    | Rattus norvegicus                                | 107            | ensembl          |
    | Saccharomy cescerevisiae(strainATCC204508/S288c) | 107            | ensembl          |

#### **Approach 2**: Run the workflow on GLDS datasets using local Ensembl reference fasta and gtf files

```bash
nextflow run NF_RCP-F_1.0.0/main.nf \
  -profile singularity \
  --gldsAccession GLDS-194 \
  --ensemblVersion 107 \
  --ref_source ensembl \
  --ref_fasta </path/to/fasta> \
  --ref_gtf </path/to/gtf>
```

> Note on ERCC References: ERCC references will still be concatenated during processing for studies that include ERCC spike-in

**Warning on Gene Annotations**: gene annotations are pulled from the current GeneLab annotation [tables](../../../Annotation_Database_Table_Generation/Gene_Annotations). This means that Ensembl/TAIR IDs that are not shared in common between the current GeneLab Ensembl version and the user supplied reference genome files  will **not** include additional gene annotations.

**Parameter Definitions:**

- `NF_RCP-F_1.0.0/main.nf` - The nextflow workflow file
- `--gldsAccession` - The GLDS accession ID to process
- `-profile` – Designates the configuration profile(s) to load
- `singularity` - A configuration profile to setup and use singularity
- `--ensemblVersion` – Indicates the Ensembl Release Version used. This is used to organize the storage of reference genome files.
- `107` - Example value, actual version should match the locally supplied gtf and fasta files.
- `--ref_source` – Indicates the reference genome source. This is used to organize the storage of reference genome and indices files.
- `ensembl` - Example reference source label.
- `--ref_fasta` – Path to a reference genome fasta file
- `--ref_gtf` – Path to a reference genome gtf file

#### **Approach 3**: Run the workflow with user-created runsheet

Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

```bash
nextflow run NF_RCP-F_1.0.0/main.nf \
  -profile singularity \
  --runsheetPath </path/to/runsheet> \
  --ensemblVersion 107 \
  --ref_source ensembl \
  --ref_fasta </path/to/fasta> \
  --ref_gtf </path/to/gtf>

```

**Parameter Definitions:**

- `NF_RCP-F_1.0.0/main.nf` - The nextflow workflow file
- `--runsheetPath` - Indicates usage of user created runsheet, replaces GeneLab accession ID
- `</path/to/runsheet>` - Path to the user created runsheet
- `-profile` – Designates the configuration profile(s) to load
- `singularity` - A configuration profile to setup and use singularity for all software
- `--ensemblVersion` – Indicates the Ensembl Release Version used. This is used to organize the storage of reference genome files.
- `107` - Example value, actual version should match the locally supplied gtf and fasta files.
- `--ref_source` – Indicates the reference genome source. This is used to organize the storage of reference genome and indices files.
- `ensembl` - Example reference source label.
- `--ref_fasta` – Path to a reference genome fasta file
- `--ref_gtf` – Path to a reference genome gtf file

**Additional Optional Arguments:**

These can be displayed as follows and include debug and storage related options not immediately useful for most users:

```bash
nextflow run NF_RCP-F_1.0.0/main.nf --help
```

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
     - VV_Logs/VV_log_VV_RAW_READS.tsv (table containing V&V flags ONLY for raw reads checks)
     - VV_Logs/VV_log_VV_TRIMMED_READS.tsv (table containing V&V flags for trimmed reads checks ONLY)
     - VV_Logs/VV_log_VV_STAR_ALIGNMENTS.tsv (table containing V&V flags for alignment file checks ONLY)
     - VV_Logs/VV_log_VV_RSEQC.tsv (table containing V&V flags for RSeQC file checks ONLY)
     - VV_Logs/VV_log_VV_RSEM_COUNTS.tsv (table containing V&V flags for RSEM raw count file checks ONLY)
     - VV_Logs/VV_log_VV_DESEQ2_ANALYSIS.tsv (table containing V&V flags for DESeq2 Analysis output checks ONLY)

Standard Nextflow resource usage logs are also produced as follows:
Further details about these logs can also found within this Nextflow documentation [page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

**Nextflow Resource Usage Logs**

   - Output:
     - Resource_Usage/execution_report_{timestamp}.html (an html report containing many useful metrics about a workflow execution including computational resources and exact workflow process commands)
     - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in your pipeline)
     - Resource_Usage/execution_trace_{timestamp}.txt (an tsv file that contains useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

---

### 6. Known Issues to Look Out For

#### Recent Nextflow versions introduced a bug that causes the workflow to fail during retrieval of files from GeneLab Repository

- [Github Issue Link](https://github.com/nextflow-io/nextflow/issues/2918)
- This will be fixed in an upcoming release of Nextflow. In the meantime, the workflow should work with Nextflow Version 21.10.6, which predates the introduction of the bug.
  - We recommend setting the environment variable 'NXF_VER=21.10.6' to allow Nextflow to automatically update/downgrade to that version on launch.
