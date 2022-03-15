#' This functions takes raw counts and labels of reference/query set to generate scGCN training input
#' @param count.list list of reference data and query data; rows are genes and columns are cells
#' @param label.list list of reference label and query label (if any), both are data frames with rownames identical with colnames of data; the first column is cell type
#' @return This function returns files saved in folders "input" & "process_data"
#' @export: all files are saved in current path
#' @examples: load count.list and label.list from folder "example_data"
#' save_processed_data(count.list,label.list)

source('data_preprocess_utility.R')

count.list <- readRDS('example_data/count.list.RDS')
label.list <- readRDS('example_data/label.list.RDS')

save_processed_data(count.list,label.list)
#四个文件，一个参考基因表达矩阵，一个待测基因表达矩阵，一个参考标签数据框，一个待测标签数据框
#> source('data_preprocess_utility.R')
#> ref_c=read.table("参考基因表达矩阵",header=T)
#> qur_c=read.table("待测基因表达矩阵",header=T)
#> count.list=list(ref_c,qur_c)
#> ref_t=read.table("ref.class.txt参考的标签数据框",header=T)
#> qur_t=read.table("qur.class.txt待测的标签数据框",header=T)
#> label.list=list(ref_t,qur_t)
#> save_processed_data(count.list,label.list)


