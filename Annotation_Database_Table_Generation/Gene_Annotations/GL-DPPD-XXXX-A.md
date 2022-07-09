# Workflow for Generating Annotation Database Tables

> **This page holds an overview and documents the generation of database gene annotation files used. Exact processing commands for specific database files that have been released are available in this repository [GLDS_Database_File_Scripts](GLDS_Database_File_Scripts) sub-directory.**  

---

**Date:** TBA
**Revision:** -  
**Document Number:** TBA

**Submitted by:**  
TBA

**Approved by:**  
TBA 

---

# Table of contents  

- [Workflow for Generating Annotation Database Tables](#workflow-for-generating-annotation-database-tables)
- [Table of contents](#table-of-contents)
- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Generate GeneLab Annotation Table File](#1-generate-genelab-annotation-table-file)

---

# Software used

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|R|4.2.0|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|`BiocManager::version()`|[https://bioconductor.org](https://bioconductor.org)| #TODO: Add Biocoductor version
|rtracklayer|1.56.0|[https://bioconductor.org/packages/release/bioc/html/rtracklayer.html](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)|
|tidyverse|1.3.1|[https://www.tidyverse.org](https://www.tidyverse.org)|
|STRINGdb|2.8.4|[https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|PANTHER.db|1.0.11|[https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html)|
|org.Hs.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)|
|org.Mm.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html)|
|org.Dm.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html)|
|org.At.tair.db|3.15.1|[https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)|

# General processing overview with example commands

> Exact processing commands for specific datasets is available in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory of this repository, as well as being provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).  

---

## 1. Generate GeneLab Annotation Table File

```bash
Rscript build-genome-annots-tab.R <MOUSE|HUMAN|ARABIDOPSIS|FLY>
```

**Parameter Definitions:**

- `<MOUSE|HUMAN|ARABIDOPSIS|FLY>` – Positional argument that specifies the organism to generate an annotation file for. Currently supports the four organisms noted.

**Input data:**

- No input files required

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations, used to add gene annotations in other workflows)
