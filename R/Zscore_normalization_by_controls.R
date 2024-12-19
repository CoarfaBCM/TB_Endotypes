# z-score normalization by controls
## exprsdf: Expression matrix with genes as row names and samples as column names.
## metadf: Annotation matrix with 1st column containing sample names, 2nd column containing "Test" or "Control" to designate the sample. The sample names should be in the same order as the column names for the exprsdf.
## Usage:
### finaldf <- zscore_norm_by_ctrls(exprsdf, metadf)
### finaldf contains the z-score normalized dataset

zscore_norm_by_ctrls <- function(exprsdf, metadf){
  
  ctrls.idx <- metadf[,2]=="Control"
  ctrls.mean <- apply(exprsdf[,ctrls.idx],1,mean)
  ctrls.sd <- apply(exprsdf[,ctrls.idx],1,sd)
  
  exprsdf <- (exprsdf-ctrls.mean)/ctrls.sd
  return(exprsdf)
}