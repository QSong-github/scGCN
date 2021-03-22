#' This functions provides the optional way to generate scGCN labels
#' @param Dat1  reference data; rows are genes and columns are cells
#' @param Dat2  query data; rows are genes and columns are cells
#' @param Lab1  label of reference data; rows are genes and columns are cells
#' @param Lab2  query data; rows are genes and columns are cells
#' @param inter_graph  inter-dataset graph in input folder

#' @return This function returns a list of two statistical scores: entropy score and enrichment score that unknown cells have higher entropy score and lower enrichment score
#' @examples: first load count.list and label.list from folder "example_data" and run 'data_preprocess.R'
#' then run 'metrics(dat1,lab1,dat2,lab2,inter_graph,cluster)'

suppressMessages(library(Seurat)); suppressMessages(library(entropy))

metrics <- function(lab1,inter_graph,clusters){
    inter.graph <- inter_graph; hc <- clusters
    
    #' ------------ check cluster enrichment in cell types ------------    
    #' Examine in clusters
    
    #' check enrichment
    eval1 <- sapply(names(table(hc)),function(name){
        ps1 <- which(hc==name)
        foo <- table(lab1[inter_graph[which(inter_graph[,2]%in%ps1),1],1])
        tem <- table(lab1)/sum(table(lab1))
        if (length(foo)>0){
            temp <- unlist(sapply(1:length(foo),function(i){
                ind <- grep(names(foo)[i],names(tem))
                r1 <- foo[i]/sum(foo)
                return (r1/tem[ind])}))
            return (max(temp)/sum(temp))
        } else {return (0)}
        })    
    
    #' check mixtureness
    eval2 <- sapply(names(table(hc)),function(name){
        ps1 <- which(hc==name)
        foo <- table(lab1[inter_graph[which(inter_graph[,2]%in%ps1),1],1])
        tem <- table(lab1)/sum(table(lab1))
        if (length(foo)>0){
            temp <- unlist(sapply(1:length(foo),function(i){
                ind <- grep(names(foo)[i],names(tem))
                r1 <- foo[i]/sum(foo)
                return (r1/tem[ind])}))
            temp2 <- entropy.empirical(temp,unit='log')
            return (temp2)
        } else {return (0) }
    })
    return (list(eval1=eval1,eval2=eval2))
}


#' This functions provides the optional way to generate scGCN labels
#' @param Dat1  reference data; rows are genes and columns are cells
#' @param Dat2  query data; rows are genes and columns are cells
#' @param Lab1  label of reference data; rows are genes and columns are cells
#' @param Lab2  query data; rows are genes and columns are cells

#' @return This function returns files saved in folders "input" & "process_data"
#' @export: all files are saved in current path
#' @examples: load count.list and label.list from folder "example_data"
#' save_processed_data(count.list,label.list)

GenerateGraph <- function(Dat1,Dat2,Lab1,K,check.unknown){
    object1 <- CreateSeuratObject(counts=Dat1,project = "1",assay = "Data1",
                                  min.cells = 0,min.features = 0,
                                  names.field = 1,names.delim = "_")
        
    object2 <- CreateSeuratObject(counts=Dat2,project = "2",assay = "Data2",
                                  min.cells = 0,min.features =0,names.field = 1,
                                  names.delim = "_")
                                  
    objects <- list(object1,object2)    
    objects1 <- lapply(objects,function(obj){
        obj <- NormalizeData(obj,verbose=F)
        obj <- FindVariableFeatures(obj,
                                       selection.method = "vst",
                                       nfeatures = 2000,verbose=F)
          obj <- ScaleData(obj,features=rownames(obj),verbose=FALSE)
          obj <- RunPCA(obj, features=rownames(obj), verbose = FALSE)
        return(obj)})
    #'  Inter-data graph  
    object.nn <- FindIntegrationAnchors(object.list = objects1,k.anchor=K,verbose=F)
    arc=object.nn@anchors
    d1.arc1=cbind(arc[arc[,4]==1,1],arc[arc[,4]==1,2],arc[arc[,4]==1,3]) 
    grp1=d1.arc1[d1.arc1[,3]>0,1:2]-1
    
    if (check.unknown){
        obj <- objects1[[2]]
        obj <- RunPCA(obj, features = VariableFeatures(object = obj),npcs=30,verbose=F)
        obj <- FindNeighbors(obj,verbose=F)
        obj <- FindClusters(obj, resolution = 0.5,verbose=F)
        hc <- Idents(obj); inter.graph=grp1+1
        scores <- metrics(lab1=Lab1,inter_graph=inter.graph,clusters=hc)
        saveRDS(scores,file='./input/statistical_scores.RDS')
    }
    #'  Intra-data graph  
    d2.list <- list(objects1[[2]],objects1[[2]])
    d2.nn <- FindIntegrationAnchors(object.list =d2.list,k.anchor=K,verbose=F)    
    d2.arc=d2.nn@anchors
    d2.arc1=cbind(d2.arc[d2.arc[,4]==1,1],d2.arc[d2.arc[,4]==1,2],d2.arc[d2.arc[,4]==1,3])
    d2.grp=d2.arc1[d2.arc1[,3]>0,1:2]-1
    final <- list(inteG=grp1,intraG=d2.grp)
    return (final)
}

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
        norm.list[[i]] <- as.matrix(Seurat:::NormalizeData.default(count.list[[i]],verbose=F))
        #' select variable features
        hvf.info <- Seurat:::FindVariableFeatures.default(count.list[[i]],selection.method='vst',verbose=F)
        hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
        hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
        var.features[[i]] <- head(rownames(hvf.info), n = 2000)
    }
    #' select variable features
    sel.features <- selectHVFeature(count.list,var.features)
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
        Seurat:::ScaleData.default(object = mat, features = hvg.features,verbose=F)})
    scale.list <- lapply(1:length(count.list),function(i){
        return (scale.list[[i]][na.omit(match(rownames(count.list[[i]]),rownames(scale.list[[i]]))),])})
    return (scale.list)}

#' This function returns a set of highly variable features
#' @param count.list list of raw counts of (1) reference data and (2) query data, rows are genes and columns are cells
#' @param var.features variable features of data list
#' @return This function returns a common set of highly variable features
#' @export
#' @examples
#' selectHVFeature(count.list,var.features)
selectHVFeature <- function(count.list,var.features,nfeatures = 2000){
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
    return (count_list_new)
}

#' This functions takes raw counts and labels of reference/query set to generate scGCN training input
#' @param count.list list of reference data and query data; rows are genes and columns are cells
#' @param label.list list of reference label and query label (if any), both are data frames with rownames identical with colnames of data; the first column is cell type
#' @return This function returns files saved in folders "input" & "process_data"
#' @export: all files are saved in current path
#' @examples: load count.list and label.list from folder "example_data"
#' save_processed_data(count.list,label.list)
save_processed_data <- function(count.list,label.list,Rgraph=TRUE,check_unknown=FALSE){
    count.list <- pre_process(count_list=count.list,label_list=label.list)
    #' save counts data to certain path: 'input'
    dir.create('input'); 
    write.csv(t(count.list[[1]]),file='input/Data1.csv',quote=F,row.names=T)
    write.csv(t(count.list[[2]]),file='input/Data2.csv',quote=F,row.names=T)

    #' optional graph: R genreated graph has minor differnce with python, user can choose the one with better performance
    if (Rgraph){
        #' use R generated graph
        new.dat1 <- count.list[[1]]; new.dat2 <- count.list[[2]]
        new.lab1 <- label.list[[1]]; new.lab2 <- label.list[[2]]
        graphs <- suppressWarnings(GenerateGraph(Dat1=new.dat1,Dat2=new.dat2,
                                                 Lab1=new.lab1,K=5,
                                                 check.unknown=check_unknown))
        write.csv(graphs[[1]],file='input/inter_graph.csv',quote=F,row.names=T)
        write.csv(graphs[[2]],file='input/intra_graph.csv',quote=F,row.names=T)
        write.csv(new.lab1,file='input/label1.csv',quote=F,row.names=F)
        write.csv(new.lab2,file='input/label2.csv',quote=F,row.names=F)        
    } else {
        #' use python generated graph
        dir.create('results')
        #' @param norm.list normalized data
        res1 <- normalize_data(count.list)
        norm.list <- res1[[1]]; hvg.features <- res1[[2]];
        #' @param scale.list scaled data
        scale.list <- scale_data(count.list,norm.list,hvg.features)
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
}
