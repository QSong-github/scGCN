# single cell graph convolutional network (scGCN)

This is a TensorFlow implementation of scGCN for leveraging differnt single cell datasets, which is described in our paper:
 
## Overview

Single-cell omics represent the fastest-growing genomics data type in the literature and the public genomics repositories. Leveraging the growing repository of labeled datasets and transferring labels from existing datasets to newly generated datasets will empower the exploration of the single-cell omics. The current label transfer methods have limited performance, largely due to the intrinsic heterogeneity among cell populations and extrinsic differences between datasets. Here, we present a robust graph artificial intelligence model, single-cell Graph Convolutional Network (scGCN), to achieve effective knowledge transfer across disparate datasets. Benchmarked with other label transfer methods on different single cell omics datasets, scGCN has consistently demonstrated superior accuracy on leveraging cells from different tissues, platforms, and species, as well as cells profiled at different molecular layers. scGCN is implemented as an integrated workflow and provided here. 

## Installation

```bash
python setup.py install
```
The general installation time is less than 10 seconds, and have been tested on mac OS and linux system. 

## Requirements
* tensorflow (>0.12)
* networkx

## Run the demo

load the example data using the data_preprocess.R script
In the example data, we include the data from Mouse (reference) and Human (query) of GSE84133 dataset. The reference dataset contains 1,841 cells and the query dataset contains more cells (N=7,264) and 12,182 genes. 
```bash
cd scGCN
Rscript data_preprocess.R # load example data 
python train.py # run scGCN
```
All output will be shown in the output_log.txt file. Performance will be shown at the bottom. 
We also provide the Seurat performance on this reference-qeury set (as in Figure 4), by run 

```
Rscript Seurat_result.R
```

## Input data

When using your own data, you have to provide 
* the raw data matrix of reference data and cell labels
* the raw data matrix of query data

## Output

The output files with scGCN predicted labels will be stored in the results folder.

## Model options 

We also provide other GCN models includidng GAT (Veličković et al., ICLR 2018), HyperGCN (Chami et al., NIPS 2019) and GWNN (Xu et al., ICLR 2019) for optional use.

## Detecting unknown cells

For the query data that have cell types not appearing in reference data, we provide a screening step in our scGCN model using two statistical metrics, entropy score and enrichment score. If certain cells in query data have higher entropy and lower enrichment, these cells should be assigned as unknown cells. Specifically, choose check_unknown=TRUE in the function 'save_processed_data' to detect unknown cells.

## Reproduction instructions

The above scripts can reproduce the quantitative results in our manuscript based on our provided data.

## Cite

Please cite our paper and the related GCN papers if you use this code in your own work:

```
Song, Q., Su, J., & Zhang, W. (2020). scGCN: a Graph Convolutional Networks Algorithm for Knowledge Transfer in Single Cell Omics. bioRxiv.
```

