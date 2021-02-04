library(Seurat)
#' @param count.list
#' @param label.list

count.list <- readRDS('./example_data/count.list.RDS')
label.list <- readRDS('./example_data/label.list.RDS')

object1 <- CreateSeuratObject(
        counts=count.list[[1]],
        project = "reference",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = label.list[[1]])
    
object2 <- CreateSeuratObject(
    counts=count.list[[2]],
    project = "query",
    assay = "RNA",
    min.cells = 0,
    min.features =0,
    names.field = 1,
    names.delim = "_",
    meta.data = label.list[[2]])
    
objs <- list(object1,object2)
objs1 <- lapply(objs,function(indrop){
    indrop <- NormalizeData(indrop)
    indrop <- FindVariableFeatures(indrop,
                                   selection.method = "vst",
                                   nfeatures = 2000)
    return(indrop)})

reference.object <- objs1[[1]]; query.object <- objs1[[2]]
    
reference.object <- ScaleData(reference.object, verbose = FALSE)
reference.object <- RunPCA(reference.object, npcs = 30, verbose = FALSE)

reference.anchors <- FindTransferAnchors(reference = reference.object, query = query.object, dims = 1:30)
predictions <- TransferData(anchorset = reference.anchors, refdata = as.factor(reference.object$type), dims = 1:30)

prediction.match <- predictions$predicted.id == query.object$type

print (sum(prediction.match)/length(prediction.match))


