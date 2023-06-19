![](https://github.com/SunXQlab/CCC-Benchmark/blob/main/Graphic-abstrct.png)  
## Introduction

We benchmark two types of CCC inference methods, one type of methods predict LR pairs based on scRNA-seq data, and another type of methods that can predict ligand/receptor-targets regulations. 

For the first benchmark, we evaluated the accuracy, stability and usability of 18 LR inference methods. **In term of accuracy**, paired ST datasets, CAGE expression/Proteomics data and sampled scRNA-seq datasets were used to benchmark the 18 methods. Firstly, 11 scRNA-seq datasets were used as input for methods to predict intercellular communication and the two defined similarity index (SI, modified Jaccard index) and rank-based similarity index (RSI) were used to compare the similarity of LR pairs predicted by methods.Furthermore, we benchmark the 18 methods using 11 paired ST datasets with the hypothesis that the values of mutual information (MI) of LR pairs are greater in the close group than that in the distant group. In addition, three PBMC datasets from 10X Genomics website were used as input for methods to predict LR pairs and CAGE expression/Proteomics data were used as pseudo gold standards to benchmark the 18 methods. **In term of stability**, we ramdomly sampled different ratios of cells in all the scRNA-seq, resulting 70 sampled datasets and 14 original datasets as input for methods. We calculated the Jaccard index of the LR pairs predicted based sampled datasets and original datasets and a stability value was defined to test the robustness of methods to sampling rates of scRNA-seq data. **In term of usability**, we recorded the running time and maximum memory usage of methods in all the 84 scRNA-seq datasets.  

For the second benchmark, 8 ST datasets were used as the input for 5 LR-Targets inference tools to predict ligand/receptor-targets regulations, and the cell line perturbation datasets were used for evaluation, involving knockout/mutant conditions for 5 receptors, and treatment conditions for 10 ligands. And the differentially expressed genes (DEGs) in each cell line perturbation dataset, were used as the ground truth of ligand/receptor-targets regulations. The score of ligand/receptor-targets predicted by different tools were compared to the differential expression status (DGEs or not DEGs) of corresponding targets to calculate AUROC and AUPRC. In addition, we also record the running time and maximum memory usage of methods in all the ST datasets. 

## Workflow
![](https://github.com/SunXQlab/CCC-Benchmark/blob/main/Workflow-figure.png)  
-   **Step0\_LRToolsFunction** contains the R/Python/Shell scripts that package the running code of 19 methods with Seurat objects as input into function.  
-   **Step1\_LRPredictionResult** contains the R/Shell scripts to run 19 methods for inferring LR pairs from the 14 scRNA-seq datasets.  
-   **Step2\_PreSTForLRBench** contains the R scripts to get the different ratios (e.g.top 10%, 20%, 30%, 40%) of cell type specific close and distant cell pairs in each dataset for the preparation of the benchmarking using mutual infomation.  
-   **Step3\_MIForLRBench** contains the R scripts to calculate MI of LR interactions predicted by methods in the different ratios of cell type specific close and distant groups and calculate DLRC index of methods in each dataset.  
-   **Step4\_SIRSIForLRBench** contains the R scripts to benchmark the similarity (SI and RSI) of the LR interactions predicted by each two methods.  
-   **Step5\_BenchBasedCAGEProteomic** contains the R scripts to benchmark the 18 LR inference methods using the CAGE expression and proteomics data.  
-   **Step6\_LRBenchSampling** contains the R/Shell scripts to run the 18 LR inference methods for inferring LR pairs from 70 sampled scRNA-seq datasets.  
-   **Step7\_LRBenchSamplingBench** contains the R/Shell scripts to calculate Jaccard index between the LR pairs predicted based on the sampled datasets and the original datasets, and record the running time and maximum memory usage of methods in each dataset.  
-   **Step8\_LRTToolsFunction** contains the R/Python/Shell scripts to run the 5 LR-Target inference methods for predicting ligand/receptor-targets using ST datasets as input.  
-   **Step9\_LRTBench** contains the R scripts to benchmark the 5 LR-Target inference methods using cell line perturbation datasets for evaluation, and record the running time and maximum memory usage of methods in each dataset.  

## Datasets
-   **scRNA-seq and ST datasets**   


<body>
    <table border="1" cellspacing="1" width="900">
        <thead>
            <tr>
                <th>Tissue (Disease)</th><th>SampleID <br>(scRNA-seq)</br></th><th>SampleID <br>(ST)</br></th><th>Literature PMID</th><th>Download URL <br>(scRNA-seq)</br></th><th>Download URL <br>(ST)</br></th><th>Evaluation purpose</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td rowspan="2">Heart Tissue (Health)</td><td>CK357</td><td>control_P7</td><td rowspan="8">35948637</td><td rowspan="8"><a href="https://zenodo.org/record/6578617">URL</a></td><td rowspan="8"><a href="https://zenodo.org/record/6580069">URL</a></td><td rowspan="2">LR interactions<br>LR-Target regulations</br></td>
            </tr>
            <tr>
                <td>CK358</td><td>control_P8</td>
            </tr>
            <tr>
                <td rowspan="3">Heart Tissue (ICM)</td><td>CK368</td><td>FZ_GT_P19</td><td rowspan="6">LR interactions</td>
            </tr>
            <tr>
                <td>CK162</td><td>FZ_GT_P4</td>
            </tr>
            <tr>
                <td>CK362</td><td>RZ_P11</td>
            </tr>
            <tr>
                <td rowspan="3">Heart Tissue (AMI)</td><td>CK361</td><td>IZ_P10</td>
            </tr>
            <tr>
                <td>CK161</td><td>IZ_P3</td>
            </tr>
            <tr>
                <td>CK165</td><td>IZ_BZ_P2</td>
            </tr>
            <tr>
                <td rowspan="2">Tumor Tissue<br>(Breast cancer)</br></td><td>CID44971</td><td>CID44971</td><td rowspan="2">34493872</td><td rowspan="2"><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078">URL</a></td><td rowspan="2"><a href="https://zenodo.org/record/4739739">URL</a></td><td rowspan="2">LR interactions<br>LR-Target regulations</br></td>
            </tr>
            <tr>
                <td>CID4465</td><td>CID4465</td>
            </tr>
            <tr>
                <td>Mouse embryo</td><td>——</td><td>Slide14</td><td>34210887</td><td>——</td><td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166692">URL</a></td><td>LR interactions</td>
            </tr>
            <tr>
                <td rowspan="3">PBMC</td><td>PBMC4K</td><td>——</td><td>——</td><td><a href="https://www.10xgenomics.com/resources/datasets/4-k-pbm-cs-from-a-healthy-donor-2-standard-2-0-1">URL</a></td><td>——</td><td rowspan="3">LR interactions</td>
            </tr>
            <tr>
                <td>PBMC6K</td><td>——</td><td>——</td><td><a href="https://www.10xgenomics.com/resources/datasets/6-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0">URL</a></td><td>——</td>
            </tr>
            <tr>
                <td>PBMC8K</td><td>——</td><td>——</td><td><a href="https://www.10xgenomics.com/resources/datasets/8-k-pbm-cs-from-a-healthy-donor-2-standard-2-0-1">URL</a></td><td>——</td>
            </tr>
            <tr>
                <td rowspan="4">Tumor Tissue<br>(Gliomas)</br></td><td>——</td><td>UKF243_T_ST</td><td rowspan="4">35700707</td><td rowspan="4">——</td><td rowspan="4"><a href="https://doi.org/10.5061/dryad.h70rxwdmj">URL</a></td><td rowspan="4">LR-Target interactions</td>
            </tr>
            <tr>
                <td>——</td><td>UKF260_T_ST</td>
            </tr>
            <tr>
                <td>——</td><td>UKF266_T_ST</td>
            </tr>
            <tr>
                <td>——</td><td>UKF334_T_ST</td>
            </tr>
        </tbody>
    </table>
</body>
  
-   **Cell line perturbation datasets**  


<body>
    <table border="1" cellspacing="1" width="900">
        <thead>
            <tr>
                <th>Datasets</th><th>Ligand/Receptor</th><th>Type</th><th>Condition</th><th>Cell Line</th><th>Disease</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td>GSE120268</td><td>AXL</td><td>receptor</td><td>Knockdown</td><td>MDA-MB-231</td><td rowspan="10">Breast Cancer</td>
            </tr>
            <tr>
                <td>GSE157680</td><td>NRP1</td><td>receptor</td><td>Knockdown</td><td>MDA-MB-231</td>
            </tr>
            <tr>
                <td rowspan="2">GSE15893</td><td>CXCR4</td><td>receptor</td><td>Mutant</td><td>MDA-MB-231</td>
            </tr>
            <tr>
                <td>CXCL12</td><td>ligand</td><td>Treatment</td><td>MDA-MB-231</td>
            </tr>
            <tr>
                <td>GSE160990</td><td>TGFB1</td><td>ligand</td><td>Treatment</td><td>MDA-MB-231</td>
            </tr>
            <tr>
                <td rowspan="3">GSE36051</td><td>DLL4(1)</td><td>ligand</td><td>Treatment</td><td>MCF7</td>
            </tr>
            <tr>
                <td>DLL4(2)</td><td>ligand</td><td>Treatment</td><td>MDA-MB-231</td>
            </tr>
            <tr>
                <td>JAG1</td><td>ligand</td><td>Treatment</td><td>MDA-MB-231</td>
            </tr>
            <tr>
                <td>GSE65398</td><td>IGF1(1)</td><td>ligand</td><td>Treatment</td><td>MCF7</td>
            </tr>
            <tr>
                <td>GSE7561</td><td>IGF1(2)</td><td>ligand</td><td>Treatment</td><td>MCF7</td>
            </tr>
            <tr>
                <td>GSE69104</td><td>CSF1R</td><td>receptor</td><td>Inhibit</td><td>TAMs</td><td rowspan="2">Gliomas</td>
            </tr>
            <tr>
                <td>GSE116414</td><td>FGFR1</td><td>receptor</td><td>Inhibit</td><td>GSLC</td>
            </tr>
            <tr>
                <td>GSE206947</td><td>EFNB2</td><td>ligand</td><td>Treatment</td><td>cardiac fibroblasts</td><td rowspan="3">Health</td>
            </tr>
            <tr>
                <td>GSE181575</td><td>TGFB1</td><td>ligand</td><td>Treatment</td><td>cardiac fibroblasts</td>
            </tr>
            <tr>
                <td>GSE123018</td><td>TGFB1</td><td>ligand</td><td>Treatment</td><td>cardiac fibroblasts</td>
            </tr>
        </tbody>
    </table>
</body>

## Tools for inferring intercellular LR pairs 

-   CellPhoneDB (Python, version: 3.0.0)
-   CellTalker (R, version: 0.0.4.9000)
-   Connectome (R, version: 1.0.1)
-   NATMI (Python)
-   ICELLNET (R, version: 1.0.1)
-   scConnect (Python, version: 1.0.3)
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
-   RNAMagnet (R, version: 0.1.0)

## Tools for predicting ligand/receptor-targets regulations

-   CytoTalk (R, version: 0.99.9)
-   NicheNet (R, version: 1.1.0)
-   stMLnet (R, version: 0.1.0)
-   MISTy (R, version: 1.3.8)
-   HoloNet (Python, version: 0.0.5)
