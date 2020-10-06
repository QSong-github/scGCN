# single cell graph convolutional network (scGCN)

This is a TensorFlow implementation of scGCN for leveraging differnt single cell datasets, which is described in our paper:
 

## Installation

```bash
python setup.py install
```

## Requirements
* tensorflow (>0.12)
* networkx

## Run the demo

load the example data using the data_preprocess.R script
In the example data, we include the data from Mouse (reference) and Human (query) of GSE84133 dataset
```bash
cd scGCN
Rscript data_preprocess.R # load example data 
python train.py # run scGCN
```
All output will be shown in the output_log.txt file

## Data

When using your own data, you have to provide 
* the raw data matrix of reference data and cell labels
* the raw data matrix of query data


## Cite

Please cite our paper if you use this code in your own work:

```
Song, Q., Su, J., & Zhang, W. (2020). scGCN: a Graph Convolutional Networks Algorithm for Knowledge Transfer in Single Cell Omics. bioRxiv.
```
