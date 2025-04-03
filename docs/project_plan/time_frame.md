## Time Frame
The table below shows the steps, tools, and estimated runtimes for each analysis.

| Date   | Hours | Process                     | Software         | Estimated Time                     | Status     |
|--------|-------|-----------------------------|------------------|------------------------------------|------------|
| 03/04  | 2     | Project Planning            | —                | 2h                                 | :white_check_mark:         |
| 07/04  | 2     | [Quality Control](#1-quality-control-and-trimming)             | FastQC           | DNA: 15 min         |            |
| 07/04  | 2     | [Short Reads preprocessing](#1-quality-control-and-trimming)    | Trimmomatic      | DNA: 5 min/sample/ (2 cores) |            |
| 07/04  | 2     | [Quality Control on Trimmed](#1-quality-control-and-trimming)   | FastQC           | DNA: 15 min       |            |
| 10/04  | 4     | [Genome Assembly](#2-genome-assembly)             | Flye             | 5h (16 cores)                      |            |
| 16/04  | 2     | [Assembly evaluation](#4-assembly-evaluation)         | QUAST and BUSCO  | 30 min                             |            |
| 20/04  | 4     | [Assembly polishing](#3-assembly-polishing)          | Pilon            | 12 h                               |            |
| 25/04  | 2     | [Assembly masking](#5-assembly-masking)            | RepeatMasker     | 1 h 40 min                         |            |
| 28/04  | 4     | [RNA mapping](#6-rna-mapping)                 | hisat2             | 15 min                             |            |
| 29/04  | 4     | [Structural annotation](#7-structural-annotation)       | BRAKER2          | 3 h 20 min                         |            |
| 02/05  | 4     | [Functional annotation](#8-functional-annotation)       | eggNOGmapper     | 2 h                                |            |
| 05/05  | 4     | [Read Counting](#9-read-counting)               | featureCounts    | 12 h                               |            |
| 07/05  | 6     | [Differential Expression](#10-differential-expression-analysis)     | DESeq2           | 10 min                             |            |
| 09/05  | 2     | [(Extra) Assembly scaffolding](#11-assembly-scaffolding-extra)| Yahs             | 1 h                                |            |
| 12/05  | 2     | [(Extra) Gene enrichment](#12-gene-enrichment-analysis-extra)     | rrvgo            | 10 min                             |            |
| 23/05  | 2     | Submit the Wiki                 | GitHub           | —                                  |            |
| 28/05  | 4     | Project Presentation        | Google Slides    | 4 h                                |            |