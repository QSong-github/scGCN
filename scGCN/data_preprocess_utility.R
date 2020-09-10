#' This function normalize count data matrix from refrence and query set
#' @param count.list list of (1) reference data and (2) query data, rows are genes and columns are cells
#' @return This function returns normalized data list 
#' @export
#' @examples
#' normalize_data(count.list)
normalize_data <- function(count.list){
    norm.list <- vector('list')
    var.features <- vector('list')
    for ( i in 1:length(count.list)){
        norm.list[[i]] <- as.matrix(Seurat:::NormalizeData.default(count.list[[i]]))
        #' select variable features
        hvf.info <- Seurat:::FindVariableFeatures.default(count.list[[i]],selection.method='vst')
        hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
        hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
        var.features[[i]] <- head(rownames(hvf.info), n = 2000)
    }
    #' select integration features
    sel.features <- selectIntegrationFeature(count.list,var.features)
    return (list(norm.list,sel.features))}

#' This function returns scaled  data matrix from refrence and query set
#' @param count.list list of raw counts of (1) reference data and (2) query data, rows are genes and columns are cells
#' @param norm.list list of normalized (1) reference data and (2) query data, rows are genes and columns are cells
#' @param hvg.features variable features of data list
#' @return This function returns scaled data list 
#' @export
#' @examples
#' scale_data(count.list,norm.list,hvg.features)
scale_data <- function(count.list,norm.list,hvg.features){
    scale.list <- lapply(norm.list,function(mat){
        Seurat:::ScaleData.default(object = mat, features = hvg.features)})
    scale.list <- lapply(1:length(count.list),function(i){
        return (scale.list[[i]][na.omit(match(rownames(count.list[[i]]),rownames(scale.list[[i]]))),])})
    return (scale.list)}

#' This function returns a set of highly variable features
#' @param count.list list of raw counts of (1) reference data and (2) query data, rows are genes and columns are cells
#' @param var.features variable features of data list
#' @return This function returns a common set of highly variable features
#' @export
#' @examples
#' selectIntegrationFeature(count.list,var.features)
selectIntegrationFeature <- function(count.list,var.features,nfeatures = 2000){
    var.features1 <- unname(unlist(var.features))
    var.features2 <- sort(table(var.features1), decreasing = TRUE)
    for (i in 1:length(count.list)) {
        var.features3 <- var.features2[names(var.features2) %in% rownames(count.list[[i]])]}    
    tie.val <- var.features3[min(nfeatures, length(var.features3))]
    features <- names(var.features3[which(var.features3 > tie.val)])
    if (length(features) > 0) {
        feature.ranks <- sapply(features, function(x) {
            ranks <- sapply(var.features, function(y) {
                if (x %in% y) {
                    return(which(x == y))
                }
                return(NULL)
            })
            median(unlist(ranks))
        })
        features <- names(sort(feature.ranks))
    }
    features.tie <- var.features3[which(var.features3 == tie.val)]
    tie.ranks <- sapply(names(features.tie), function(x) {
        ranks <- sapply(var.features, function(y) {
            if (x %in% y) {return(which(x == y))}
            return(NULL)
        })
        median(unlist(ranks))
    })
    features <- c(features, names(head(sort(tie.ranks), nfeatures - length(features))))
    return(features)
}

#' This functions takes refrence data and labels to identify variable gene features
#' @param data reference data; rows are genes and columns are cells
#' @param label data frame with rownames identical with colnames of data; the first column is cell type
#' @param nf number of variable features
#' @return This function returns variable gene features using ANOVA
#' @export
#' @examples
#' select_feature(data=reference.data,label=reference.label)
select_feature <- function(data,label,nf=2000){
    M <- nrow(data); new.label <- label[,1]
    pv1 <- sapply(1:M, function(i){
        mydataframe <- data.frame(y=as.numeric(data[i,]), ig=new.label)
        fit <- aov(y ~ ig, data=mydataframe)
        summary(fit)[[1]][["Pr(>F)"]][1]})
    names(pv1) <- rownames(data)
    pv1.sig <- names(pv1)[order(pv1)[1:nf]]
    egen <- unique(pv1.sig)
    return (egen)
}

#' This functions takes refrence data and labels to identify variable gene features
#' @param count_list list of reference data and query data; rows are genes and columns are cells
#' @param label_list list of reference label and query label (if any), both are data frames with rownames identical with colnames of data; the first column is cell type
#' @return This function returns four elements, including the normalized data, scaled data, hvg features, and selected features
#' @export
#' @examples
#' pre_process(count_list,label_list)
pre_process <- function(count_list,label_list){
    sel.features <- select_feature(count_list[[1]],label_list[[1]])
    count_list_new <- list(count_list[[1]][sel.features,],count_list[[2]][sel.features,])
    res1 <- normalize_data(count_list_new)
    norm_list <- res1[[1]]; hvg_features <- res1[[2]]; 
    #' scale features (use normalized data)
    #' @param scale.list scaled data
    scale_list <- scale_data(count_list_new,norm_list,hvg_features)
    return (list(count_list_new,norm_list,scale_list,hvg_features))
}

#' This functions takes raw counts and labels of reference/query set to generate scGCN training input
#' @param count.list list of reference data and query data; rows are genes and columns are cells
#' @param label.list list of reference label and query label (if any), both are data frames with rownames identical with colnames of data; the first column is cell type
#' @return This function returns files saved in folders "input" & "process_data"
#' @export: all files are saved in current path
#' @examples: load count.list and label.list from folder "example_data"
#' save_processed_data(count.list,label.list)
save_processed_data <- function(count.list,label.list){
    res1 <- pre_process(count_list=count.list,label_list=label.list)

    count.list <- res1[[1]];
    norm.list <- res1[[2]];
    scale.list <- res1[[3]];
    hvg.features <- res1[[4]]

    #' save counts data to certain path: 'input'
    dir.create('input'); dir.create('results')
    write.csv(t(count.list[[1]]),file='input/Data1.csv',quote=F,row.names=T)
    write.csv(t(count.list[[2]]),file='input/Data2.csv',quote=F,row.names=T)
    #' save processed data to certain path: 'process_data'
    outputdir <- 'process_data'; dir.create(outputdir)
    write.csv(hvg.features,file=paste0(outputdir,'/sel_features.csv'),quote=F,row.names=F)

    N <- length(count.list)
    for (i in 1:N){
        df = count.list[[i]]
        if (!dir.exists(paste0(outputdir,'/count_data'))){dir.create(paste0(outputdir,'/count_data'))}
        file.name=paste0(outputdir,'/count_data/count_data_',i,'.csv')
        write.csv(df,file=file.name,quote=F)
    }

    for (i in 1:N){
        df = label.list[[i]]
        if (!dir.exists(paste0(outputdir,'/label_data'))){dir.create(paste0(outputdir,'/label_data'))}
        file.name=paste0(outputdir,'/label_data/label_data_',i,'.csv')
        write.csv(df,file=file.name,quote=F)
    }
    
    for (i in 1:N){
        df = norm.list[[i]]
        if (!dir.exists(paste0(outputdir,'/norm_data'))){dir.create(paste0(outputdir,'/norm_data'))}
        file.name=paste0(outputdir,'/norm_data/norm_data_',i,'.csv')
        write.csv(df,file=file.name,quote=F)
    }
    
    for (i in 1:N){
        df = scale.list[[i]]
        if (!dir.exists(paste0(outputdir,'/scale_data'))){dir.create(paste0(outputdir,'/scale_data'))}
        file.name=paste0(outputdir,'/scale_data/scale_data_',i,'.csv')
        write.csv(df,file=file.name,quote=F)
    }
}





