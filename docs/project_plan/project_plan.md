## Goals

- Assemble one chromosome of the *Niphotrichum japonicum* genome using Nanopore long reads and Illumina short reads.  
- Functionally annotate the assembled chromosome.  
- Identify differentially expressed genes (DEGs) under heat stress.  
- *(Extra)* Incorporate Hi-C data to improve assembly scaffolding.  
- *(Extra)* Perform gene enrichment analysis for unique and/or differentially expressed genes.


## Type of sample

For genome assembly, the sample consists of gametophytes of *Niphotrichum japonicum*. For transcriptome analysis, gametophyte samples were collected under both control conditions and after exposure to heat stress.

## Type of data

- Nanopore long reads for assembly of chromosome 3
- Illumina short reads for polishing the assembled sequence
- Illumina RNA-seq reads for transcriptome analysis under heat stress and - control conditions
- (Extra) Illumina Hi-C reads for scaffolding

## Analyses and Software

The diagram below shows the main steps of the project, organized into genome assembly (green), transcriptome and expression analysis (blue), and extra analyses (red).

![Workflow of the project](/docs/images/workflow.jpeg)

The steps below are listed in the order they will be carried out during the project.

### 1. Quality control and trimming
Illumina DNA reads will be checked for quality with FastQC before and after trimming. Trimmomatic will be used to remove adapters and low-quality bases.

### 2. Genome assembly
Nanopore long reads will be assembled into contigs using Flye.

### 3. Assembly polishing
Illumina DNA reads will be used to polish the assembly and correct sequencing errors using Pilon.

### 4. Assembly evaluation
The quality of the assembly will be assessed using QUAST and BUSCO.

### 5. Assembly masking
RepeatMasker will be used to identify and mask repetitive regions in the assembly.

### 6. RNA mapping
Trimmed RNA-seq reads will be aligned to the genome using STAR.

### 7. Structural annotation
BRAKER2 will be used to predict gene structures using the masked genome and RNA-seq read alignments.

### 8. Functional annotation
Predicted genes will be annotated with eggNOGmapper to assign putative functions.

### 9. Read counting
Trimmed RNA-seq reads will be mapped to the genome using STAR, and gene counts will be generated with featureCounts.

### 10. Differential expression analysis
Gene counts will be analyzed with DESeq2 to identify genes that are differentially expressed between control and heat stress conditions.

### 11. Assembly scaffolding (Extra)
Hi-C reads will be used with Yahs to scaffold the assembly and improve its structure.

### 12. Gene enrichment analysis (Extra)
rrvgo will be used to perform gene enrichment analysis on unique or differentially expressed genes to explore their functional roles.

## Time Frame
The table below shows the steps, tools, and estimated runtimes for each analysis.

| Date   | Hours | Process                     | Software         | Estimated Time                     | Status     |
|--------|-------|-----------------------------|------------------|------------------------------------|------------|
| 03/04  | 2     | Project Planning            | —                | 2h                                 | :white_check_mark:         |
| 07/04  | 2     | [Quality Control](#1-quality-control-and-trimming)             | FastQC           | DNA: 15 min         |            |
| 08/04  | 2     | [Short Reads preprocessing](#1-quality-control-and-trimming)    | Trimmomatic      | DNA: 5 min/sample/ (2 cores) |            |
| 09/04  | 2     | [Quality Control on Trimmed](#1-quality-control-and-trimming)   | FastQC           | DNA: 15 min       |            |
| 10/04  | 4     | [Genome Assembly](#2-genome-assembly)             | Flye             | 5h (16 cores)                      |            |
| 22/04  | 4     | [Assembly polishing](#3-assembly-polishing)          | Pilon            | 12 h                               |            |
| 22/04  | 2     | [Assembly evaluation](#4-assembly-evaluation)         | QUAST and BUSCO  | 30 min                             |            |
| 25/04  | 2     | [Assembly masking](#5-assembly-masking)            | RepeatMasker     | 1 h 40 min                         |            |
| 28/04  | 4     | [RNA mapping](#6-rna-mapping)                 | hisat2             | 15 min                             |            |
| 29/04  | 4     | [Structural annotation](#7-structural-annotation)       | BRAKER2          | 3 h 20 min                         |            |
| 02/05  | 4     | [Functional annotation](#8-functional-annotation)       | eggNOGmapper     | 2 h                                |            |
| 14/05  | 4     | [Read Counting](#9-read-counting)               | featureCounts    | 12 h                               |            |
| 20/05  | 6     | [Differential Expression](#10-differential-expression-analysis)     | DESeq2           | 10 min                             |            |
| 25/05  | 4     | Prepare Presentation        | Google Slides    | 4 h                                |            |
| 25/05  | 2     | Refine Wiki                 | GitHub           | —                                  |            |
| 28/05  | 2     | [(Extra) Assembly scaffolding](#11-assembly-scaffolding-extra)| Yahs             | 1 h                                |            |
| 28/05  | 2     | [(Extra) Gene enrichment](#12-gene-enrichment-analysis-extra)     | rrvgo            | 10 min                             |            |

## Data Management
The project is organized into separate folders for data, analyses, code, and documentation. Raw data is stored in the `data/` folder using symbolic links, each analysis step has its own subfolder in `analyses/`, and scripts are stored in `code/`.

![repository structure](/docs/images/folder_structure.png)

## References
- Zhou, X., Peng, T., Zeng, Y., Cai, Y., Zuo, Q., Zhang, L., Dong, S., & Liu, Y. 2023. Chromosome-level genome assembly of *Niphotrichum japonicum* provides new insights into heat stress responses in mosses. *Frontiers in Plant Science* 14: 1271357; [https://doi.org/10.3389/fpls.2023.1271357](https://doi.org/10.3389/fpls.2023.1271357)

