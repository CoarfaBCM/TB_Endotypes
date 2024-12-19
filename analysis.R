# box link to download "input_data" folder containing input data files used for the original analysis:
## https://bcm.box.com/s/5lxqk347vvwhkmeh34ml7jybejrpmh9y

# box link to download "ClusteringTree" folder containing clustering tree results from the original analysis:
## https://bcm.box.com/s/y0accowv3zm3j8a72nokwfpkcsc8ty6v

# reading in discovery cohort metadata
metadf <- openxlsx::read.xlsx("./input_data/Metadata-AllBalanced-MA.xlsx")

# reading in discovery cohort combat normalized expression data
exprsdf_cbat_ma <- readRDS("./input_data/exprs_cbat-AllBalanced-MA.rds")

# creating a new metadata dataframe with sample ID, test or control group, and batch number
metadata <- data.frame("ID"=metadf$GSM.code[match(colnames(exprsdf_cbat_ma),metadf$GSM.code)],
                       "Group"=metadf$Group[match(colnames(exprsdf_cbat_ma),metadf$GSM.code)],
                       "Set"=metadf$Batch[match(colnames(exprsdf_cbat_ma),metadf$GSM.code)],
                       stringsAsFactors = F)

# adding GSEs as batch names
metadata$Set <- unname(sapply(metadata$ID, function(x){metadf$GSE[metadf$GSM.code==x]}))

# reorder samples to match sample order between columns of expression data and rows of metadata
exprsdf_cbat_ma <- exprsdf_cbat_ma[, metadata$ID]
identical(colnames(exprsdf_cbat_ma), metadata$ID)

# clustering Tree
source("R/ClusteringTree.R")

# generating the clustering tree
ClusteringTree(outdir = "./ClusteringTree/", exprs_df = exprsdf_cbat_ma, meta_df = metadata, res.seq = seq(0,1,0.2), initrun = T)

# loading in clustering information
my_obj <- readRDS("./ClusteringTree/FullSeuratObj_0-1-0.2.rds")
clusteringdf <- my_obj@meta.data

# loading UQ normalized expression data from the Borstel cohort (validation cohort)
exprsdf_uq_bor <- readRDS("./input_data/Exprs_UQ-combined.rds")

# loading UQ normalized RNA-Seq expression data (validation cohort)
exprsdf_uq_rna <- readRDS("./input_data/exprs_UQ-AllBalanced-RNASeq.rds")

# reducing discovery cohort to genes common between discovery, borstel and rna-seq cohorts
exprsdf_training <- exprsdf_cbat_ma[rownames(exprsdf_cbat_ma) %in% rownames(exprsdf_uq_bor), ]
exprsdf_training <- exprsdf_training[rownames(exprsdf_training) %in% rownames(exprsdf_uq_rna), ]
identical(colnames(exprsdf_training), metadata$ID)

# z-score normalization
source("./R/Zscore_normalization_by_controls.R")
exprsdf_training <- zscore_norm_by_ctrls(exprsdf = exprsdf_training, metadf = metadata)

# using cluster assignment for samples at resolution 0.4
cluster.list <- clusteringdf[,"RNA_snn_res.0.4"]
metadata$Cluster <- sapply(metadata$ID,
                           function(x){
                             if(x %in% rownames(clusteringdf)){as.numeric(cluster.list[rownames(clusteringdf)==x])-1}else{99}
                             })
exprsdf_training_tb <- data.frame(Cluster=as.factor(metadata$Cluster[metadata$Cluster!=99]), t(exprsdf_training[, metadata$Cluster!=99]))

# RF training

# create classification output directory if it does not exist
outdir <- "./Classification_0.4/"
if (!dir.exists(outdir)){dir.create(outdir, recursive = TRUE)}

library(randomForest)
# myRF <- randomForest(Cluster~., data = exprsdf_training_tb, importance = T)
# 
# myimp <- data.frame(Gene=row.names(myRF$importance), myRF$importance, stringsAsFactors = F)
# myimp <- myimp[order(myimp$MeanDecreaseGini, decreasing = T), ]
# write.xlsx(myimp[,-c(2,3)], "./Classification_0.4/RF_Importancedf_AllBalanced_0.4_commonDBV.xlsx", rownames=F)
myimp <- read.xlsx("./Classification_0.4/RF_Importancedf_AllBalanced_0.4_commonDBV.xlsx")

# # building RF models with a range of number of genes and saving their OOB error rate estimate
# myseq <- c(seq(10,100,10), 200, 500, seq(1000,ncol(exprsdf_training_tb)-1,1000), ncol(exprsdf_training_tb)-1)
# oob.err <- c()
# allRFmodels <- list()
# 
# for (numgenes in myseq) {
#   myRF <- randomForest(Cluster ~ ., data = exprsdf_training_tb[,c("Cluster", myimp$Gene[1:numgenes])], importance = T)
#   oob.err <- c(oob.err,mean(predict(myRF)!= exprsdf_training_tb$Cluster))
#   allRFmodels[[numgenes]] <- myRF
# }
# ooberrdf <- data.frame("Genes"=myseq,"OOBErrorRate"=paste0(round(oob.err*100, 2),"%"), row.names = NULL)
# ooberrdf$Rank <- rank(oob.err, ties.method = "min")
# write.xlsx(ooberrdf,"./Classification_0.4/OOBerrorrate.xlsx", overwrite = T)

# # after investigating OOBerrorrate.xlsx, we pick the 50 gene model as the best one to move forward with
# # running validation on rna-seq and borstel data using this 50 gene RF classifier model
# numgenes <- 50
# 
# myRF <- allRFmodels[[numgenes]]
# saveRDS(myRF, paste0("./Classification_0.4/RF_Classifier_",numgenes,"genes_AllBalanced_0.4_commonDBV.rds"))


# loading borstel metadata
metadata_bor <- read.xlsx("./input_data/meta_Borstel_TB_HC_Combined-Batch.xlsx")
metadata_bor$Set <- metadata_bor$Set %>% gsub(1, "GIC",.) %>% gsub(2, "GVC",.) %>% gsub(3, "RVC",.) %>% as.factor()

identical(colnames(exprsdf_uq_bor), metadata_bor$ID)

# loading rna-seq metadata
metadf_rna <- read.xlsx("./input_data/Metadata-AllBalanced-RNASeq.xlsx")
metadata_rna <- data.frame("ID"=metadf_rna$GSM.code[match(colnames(exprsdf_uq_rna),metadf_rna$GSM.code)],
                           "Group"=metadf_rna$Group[match(colnames(exprsdf_uq_rna),metadf_rna$GSM.code)],
                           "Set"=metadf_rna$Batch[match(colnames(exprsdf_uq_rna),metadf_rna$GSM.code)])
metadata_rna$Set <- as.factor(metadata_rna$Set)
identical(colnames(exprsdf_uq_rna), metadata_rna$ID)

# validation - borstel
# Loading RF classifier model
myRF <- readRDS(paste0("./Classification_0.4/RF_Classifier_",numgenes,"genes_AllBalanced_0.4_commonDBV.rds"))

myRFpred <- function(exprs_df, meta_df, cohort, numgenes) {
  source("./R/Zscore_normalization_by_controls.R")
  
  topgenes <- myimp$Gene[1:numgenes]
  # Z-Score Normalization
  # range(exprsdf)
  testingdf <- zscore_norm_by_ctrls(exprs_df, meta_df)
  # range(exprsdf)

  ## Restricting testingdf to top 50 genes and picking out only Test samples
  ## Check the conversion of "-" to "." in gene names to keep uniformity
  testingdf <- t(testingdf[topgenes,meta_df$Group=="Test"])

  # Testing
  myprob_RF <- predict(myRF, testingdf, type = "prob")
  # table(round(myprob_RF[,2]))
  resultsdf <- data.frame(Sample=row.names(myprob_RF), Cluster=round(myprob_RF[,2]))
  meta_df$Cluster <- sapply(meta_df$ID,function(x){ifelse(x %in% resultsdf$Sample, resultsdf$Cluster[resultsdf$Sample==x], 99)})
  table(meta_df$Group,meta_df$Cluster)
  
  write.table(meta_df, paste0("./Classification_0.4/RF_FullClusterPred_",numgenes,"genes_",cohort,"_AllBalanced_0.4_commonDBV.xls"), row.names = F, sep = "\t")
}
myRFpred(exprs_df = exprsdf_uq_bor, meta_df = metadata_bor, cohort = "Borstel", numgenes = 50)
myRFpred(exprs_df = exprsdf_uq_rna, meta_df = metadata_rna, cohort = "RNA", numgenes = 50)
