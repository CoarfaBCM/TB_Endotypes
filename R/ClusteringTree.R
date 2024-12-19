# Clustering Tree Analysis

ClusteringTree <- function(outdir=getwd(),
                                 exprs_df,
                                 meta_df,
                                 res.seq,
                                 initrun=T,
                                 PCAdims=50,
                                 prefix=NULL,
                                 suffix=NULL,
                                 fc=c(1, 1.5, 2),
                                 batch=F,
                                 limma=F,
                                 perplexity = 30) {
  
  # Code to install bioconductor packages
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install(c("limma"))
  
  # Loading required packages and installing ones not present
  list.of.packages <- c("Seurat","RColorBrewer","openxlsx","ggfortify","clustree","tidyverse","ggplot2")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  # Making sure exprsdf and metadf samples are in the same order
  colnames(meta_df) <- c("ID", "Group", "Set")[1:ncol(meta_df)]
  if (isTRUE(batch)) {meta_df$Set <- as.factor(meta_df$Set)}
  meta_df <- meta_df[match(colnames(exprs_df), meta_df$ID),]
  
  #Taking just the Test samples
  newexprsdf <- exprs_df[,meta_df$ID[meta_df$Group=="Test"]]
  
  ##Seurat
  my_obj <- CreateSeuratObject(newexprsdf)
  
  my_obj <- SetAssayData(my_obj, layer = "data", assay = "RNA", new.data = newexprsdf)
  
  my_obj <- FindVariableFeatures(my_obj)
  
  my_obj <- ScaleData(my_obj)
  
  my_obj <- RunPCA(my_obj)
  
  # create output directory if it does not exist
  if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)}
  
  pdf(paste0(outdir,"/ElbowPlot_",PCAdims,"dims.pdf"),height = 10, width = 10)
  print(ElbowPlot(my_obj, ndims = PCAdims))
  dev.off()
  
  my_obj <- FindNeighbors(my_obj,
                          reduction = "pca",
                          dims = 1:PCAdims,
                          verbose = 1)
  
  mylist <- list()
  samplesbycluster <- list()
  
  # Clustering Tree Loop
  for (res.ind in 1:length(res.seq)) {
    # Clustering
    cat("\n","######## Res = ", res.seq[res.ind], " ########","\n")
    my_obj <- FindClusters(my_obj, resolution = res.seq[res.ind], verbose = 1)
    
    # Saving sample clustering and number of samples per cluster at each resolution
    if(batch){
      mylist[[res.ind]] <- data.frame("Sample"=meta_df$ID,
                                      "Cluster"=sapply(meta_df$ID, function(x){if (x %in% rownames(my_obj@meta.data)) {as.numeric(my_obj@meta.data[rownames(my_obj@meta.data) == x, "seurat_clusters"])-1} else {99}}),
                                      "Batch"=meta_df$Set)
    } else {
      mylist[[res.ind]] <- data.frame("Sample"=meta_df$ID,
                                      "Cluster"=sapply(meta_df$ID, function(x){if (x %in% rownames(my_obj@meta.data)) {as.numeric(my_obj@meta.data[rownames(my_obj@meta.data) == x, "seurat_clusters"])-1} else {99}})) 
    }
    samplesbycluster[[res.ind]] <- as.data.frame(table(as.numeric(my_obj@meta.data[,"seurat_clusters"])-1))
    colnames(samplesbycluster[[res.ind]]) <- c("Cluster",paste0("Res_",res.seq[res.ind]))
  }
  samplesbycluster <- samplesbycluster %>% reduce(full_join, by="Cluster")
  names(mylist) = paste0("Res",res.seq)
  
  write.xlsx(mylist, paste0(outdir,"/SampleClusterNum_AllRes.xlsx"), row.names=F, overwrite = T)
  write.xlsx(samplesbycluster, paste0(outdir,"/NumSamplesbyCluster_AllRes.xlsx"), col.names = T, row.names=F, overwrite = T)
  
  saveRDS(my_obj, paste0(outdir,"/FullSeuratObj_",min(res.seq),"-",max(res.seq),"-",res.seq[2]-res.seq[1],".rds"))
  
  #png("clustering_tree.png", width = 8, height = 8, units = "in", pointsize = 14, bg = "white", res =300)
  pdf(paste0(outdir,"/clusteringtree_",min(res.seq),"-",max(res.seq),"-",res.seq[2]-res.seq[1],".pdf"),height = 10, width = 10)
  print(clustree(my_obj@meta.data, prefix = "RNA_snn_res.", node_label = "size") + labs(color="Res") + guides(fill="none"))
  dev.off()
  
  if (initrun) {
    cat("\n#### Initial run is complete and clustering tree is available. ####\n")
    return(invisible(NULL))
  }
}
