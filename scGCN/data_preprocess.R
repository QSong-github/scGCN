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

