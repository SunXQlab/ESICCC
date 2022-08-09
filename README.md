![](https://github.com/SunXQlab/CCC-Benchmark/blob/main/workflow.png)

## Introduction

We benchmark two types of CCC inference methods/tools, one type of methods predict LR pairs based on scRNA-seq
data, and another type of methods that can predict ligand/receiver-targets regulations. 

For the first benchmark, scRNA-seq data is used as the input data for
different tools to predict intercellular communication. We define a
similarity index (modified Jaccard index) to compare the similarity of
LR pairs predicted by tools. In addition, we hypothesize that close cell
pairs are more likely to communicate than distant cell pairs. Therefore,
we assume that the values of mutual information (MI) and Pearson correlation
coefficient (PCC) of LR pairs are greater in the close group than that in the
distant group. Using the coordinate information in the ST data, we
divide the cell pairs into close and distant groups, and then calculate
the MI and PCC values of LR pairs predicted by different tools in the close and distant
groups.

For the second benchmark, ST data is used as the input data for
different tools to predict ligand/receiver-targets regulations, and the cell line
perturbation datasets are used for evaluation, involving knockout/mutant
conditions for 3 receptors, and treatment conditions for 5 ligands. And
the differentially expressed genes (DEGs) in each cell line perturbation dataset, are used as the ground truth of
ligand/receptor-targets regulations. The score of
ligand/receptor-targets predicted by different tools were compared to
the differential expression status (DGEs or not DEGs) of corresponding
targets to calculate AUCROC and AUCPR.

## Workflow

1.  **step1\_data\_process** contains the R script to perform data
    analysis of scRNA-seq and ST data.
2.  **step2\_sc\_tools\_benchmark** contains the R script to run 19
    methods/tools for inferring LR pairs from the scRNA-seq dataset.
3.  **step3\_run\_lr\_benchmark** contains the R script to benchmark the
    above 18 methods/tools (except for cell2cell) using mutual
    information, Pearson correlation coefficient and the similarity index.
4.  **step4\_tools\_intra\_benchmark** contains the R script to run and
    benchmark the 5 methods/tools for predicting ligand/receptor-targets
    using ST dataset for input and cell line perturbation datasets for
    evaluation.

## Datasets

-   **scRNA-seq dataset** can be downloaded from the Gene Expression
    Omnibus (GEO) repository with the accession number GSE176078.
-   **ST dataset** can be downloaded from Zenodo data repository
    (<https://doi.org/10.5281/zenodo.4739739>).
-   **Cell line perturbation datasets** can be download from GEO
    repository with the accession number GSE15893, GSE120268, GSE157680,
    GSE7561, GSE15893, GSE36051, GSE65398 and GSE160990.

## Tools for inferring intercellular LR pairs 

-   CellPhoneDB (Python, version: 3.0.0)
-   CellTalker (R, version: 0.0.4.9000)
-   Connectome (R, version: 1.0.1)
-   NATMI (Python)
-   ICELLNET (R, version: 1.0.1)
-   scConnect (Python, version: 1.0.3)
-   Cellinker (Webserver)
-   CellChat (R, version: 1.4.0)
-   SingleCellSignalR (R, version: 1.2.0)
-   CytoTalk (R, version: 0.99.9)
-   CellCall (R, version: 0.0.0.9000)
-   scSeqComm (R, version: 1.0.0)
-   NicheNet (R, version: 1.1.0)
-   Domino (R, version: 0.1.1)
-   scMLnet (R, version: 0.2.0)
-   PyMINEr (Python, version: 0.10.0)
-   iTALK (R, version: 0.1.0)
-   cell2cell (Python, version: 0.5.10)

## Tools for predicting ligand/receptor-targets regulations

-   CytoTalk (R, version: 0.99.9)
-   NicheNet (R, version: 1.1.0)
-   stMLnet (R, version: 0.1.0)
-   MISTy (R, version: 1.3.8)
-   HoloNet (Python, version: 0.0.5)
