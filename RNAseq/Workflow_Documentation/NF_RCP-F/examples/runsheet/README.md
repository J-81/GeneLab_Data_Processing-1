# RunSheet Specification

## Description

1. The RunSheet is a csv file that contains the metadata required for processing bulkRNASeq datasets through GeneLab's processing workflow.

## Examples

1. [Runsheet for GLDS-194](paired_end_runsheet/GLDS-194_bulkRNASeq_v1.csv) (paired end dataset)
2. [Runsheet for GLDS-48](single_end_runsheet/GLDS-48_bulkRNASeq_v1.csv) (single end dataset)

## Specific File Usage

1. Used in the Nextflow Workflow to describe metadata required for processing
2. Used in the differential gene expression to describe experimental groups

## Required columns

| Column Name                            | Type                      | Description                                                                                                                                            | Example                       |
|----------------------------------------|---------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------|
| Sample Name                            | string                    | Sample Name used during the processing.                                                                                                                | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |
| has_ERCC                               | bool                      | Set to True if ERCC spike-ins are included. This ensures ERCC normalized DGE is performed in addition to normal DGE                                    | True                          |
| paired_end                             | bool                      | Set to True if the library layout is paired end.                                                                                                       | False                         |
| organism                               | string                    | Organism name used to map to the appropriate gene annotations file. Supported organisms can be found in organisms.csv.                                 | Mus musculus                  |
| read1_path                             | string(url or local path) | Location of the raw reads file. In paired end layouts, this denotes the forward reads file.                                                            | /my/data/sample_1.fastq.gz    |
| read2_path                             | string(url or local path) | Location of the raw reads file. In paired end layouts, this denotes the reverse reads file. For single end layouts, this column should be omitted.                                                            | /my/data/sample_2.fastq.gz    |
| Factor Value[<name, e.g. Spaceflight>] | string                    | A set of columns denoting the experimental group for the sample.  In the simplest form, a column named 'Factor Value[group]' is sufficient.            | Space Flight                  |
| Original Sample Name                   | string                    | Used to map processing sample name to original sample name. This is often identical except in cases where the original name includes space characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |