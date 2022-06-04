scRNA_GCN流程是在scGCN流程的基础上进行的升级和完善，支持单细胞RNA基因表达矩阵作为输入，用以跨物种之间的细胞类型比较。
## 概述

基于卷积神经网络（Convolutional Neural Networks）的跨物种单细胞卷积神经网络细胞类型注释软件是一个linux平台的python软件，旨在帮助没有过多的编程知识的用户在下游分析单细胞转录组测序（scRNA-seq）及相关组学数据。软件提供了一种方便的算法来对不同物种的基因表达相似性进行分析，从而在无需知道跨物种之间的物种特异性marker的情况下，根据参考物种的细胞类型注释对待测物种的细胞类型进行注释。并且，在细胞类型注释完成后，还能简单地比较分析不同物种细胞类型的基因表达差异。软件操作难度低，利于没有编程基础的研究人员学习并使用。

## 需要安装的包
* setuptools >= 40.6.3
* numpy >= 1.15.4
* tensorflow >= 1.15.0
* networkx >= 2.2
* scipy >= 1.1.0

## 流程的安装

下载流程：
```
git clone https://github.com/Dee-chen/scGCN/tree/cdy
```
安装依赖:

```bash
python setup.py install
```


## 运行示例

在示例文档中，有两个物种的相关基因表达矩阵数据
```bash
cd scGCN
Rscript data_preprocess.R # 加载和预处理示例数据
python train.py # 运行主程序
```
```
Rscript Seurat_result.R
```

## 输入数据

当使用自定义数据时，需要提供一个待测物种的基因表达矩阵，和一个参考物种的基因表达矩阵以及对应的细胞类型标签数据
* 基因表达矩阵是数据框结构，行是细胞，列是基因
* 细胞类型标签是单列数据框结构，以"type"作为单列数据的表头

## 输出

输出的数据信息会在results文件夹中，可视化结果则在当前运行的主程序目录下

## 差异分析
在原始scGCN的基础上，增加了初步的差异分析功能，用以可视化跨物种同一细胞类型的主要差异基因

## 引用

原始scGCN的文章引用信息：
```
Song, Q., Su, J., & Zhang, W. (2020). scGCN: a Graph Convolutional Networks Algorithm for Knowledge Transfer in Single Cell Omics. bioRxiv.
```

