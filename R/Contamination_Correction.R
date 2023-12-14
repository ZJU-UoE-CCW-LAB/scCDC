#' A simple ROC function to obtain the AUROC value
#'
#' Calculate the AUROC based on the values representing cellular expression levels and the
#' corresponding cell type (two types) information.
#'
#' @param exp the values representing cellular expression levels of a certain gene
#' @param class the cell types of the corresponding cells
#'
#' @return an AUROC (Area Under the Receiver Operating Characteristics) value
#'
#' @import pROC
#'
simple_roc <- function(exp,class){
  # match and sort
  names(exp) <- class
  exp_class <- sort(exp)
  # draw the roc
  curve <- roc(response=names(exp_class), predictor = exp_class,direction='<')
  AUC <- curve[['auc']]
  return(AUC)
}

#' AUROC calculation between the cluster with the lowest expression and each other cluster
#'
#' Calculate the AUROC values based on a SeuratObject and a selected gene. The cell cluster 
#' with the lowest expression level will be compared with each other cell cluster (including 
#' itself) and then generate the corresponding AUROC value. These values will output with 
#' the cluster names.
#'
#' @param object a clustered SeuratObject
#' @param gene a gene within the input SeuratObject
#' @param qualified_cls a vector of qualified cluster names
#'
#' @return a named vector of the AUROC values
#'
#' @import Seurat pROC
#'
Cal_AUCs <- function(object,gene,qualified_cls){
  # obtain the average expression level of each cluster
  ave_exp <- AverageExpression(object, features = gene)
  exps <- ave_exp[[1]]
  # sort the cluster with the average expression level
  exp_sort <- as.vector(exps)
  names(exp_sort) <- colnames(exps)
  exp_sort <- sort(exp_sort)
  cls <- names(exp_sort)
  # find the cluster with the lowest average expression level 
  # (the first eGCG_negative cluster)
  eGCG_neg <- intersect(cls,qualified_cls)[1]
  # calculate the auc between each cluster with the first eGCG_negative cluster
  aucs <- c()
  for (cluster in cls) {
    
    cls_cells <- WhichCells(object, idents = cluster)
    neg_cells <- WhichCells(object, idents = eGCG_neg)
    
    counts_exp <- GetAssayData(object, slot = 'counts')
    counts_exps <- counts_exp[gene, c(cls_cells, neg_cells)]
    
    class <- c(rep(1, length(cls_cells)), rep(0, length(neg_cells)))
    
    # apply the simple roc function to obtain the auc value
    auc_val <- simple_roc(counts_exps, class)
    aucs <- c(aucs, auc_val)
  }
  
  names(aucs) <- cls
  ###
  return(list(aucs,eGCG_neg))
}

#' Find the threshold (best cutting-point) using the Youden index
#'
#' Calculate the threshold (best cutting-point) based on the Youden indexes between the
#' count values of the eGCG_postive cluster with the lowest expression and that of all 
#' eGCG_negative cluster. The best cutting-point will be selected when Younden index is
#' maximized. The eGCG_positive and eGCG_negative clusters will be identified based on
#' the input eGCG_aucs and the user-specified auc_thres. The cluster with the auc value
#' smaller than the auc_thres will be identified as eGCG_negative cluster.
#'
#' @param object a clustered SeuratObject
#' @param gene a gene within the input SeuratObject
#' @param eGCG_aucs a named vector of the AUROC values between the each cluster and the cluster with the lowest expression level
#' @param auc_thres the AUROC threshold to determine the boundary between eGCG_positive and eGCG_negative clusters
#'
#' @return the count value of the threshold (best cutting-point)
#'
#' @import Seurat pROC ddpcr
#'
Cal_thres <- function(object,gene,eGCG_aucs,auc_thres){
  # setting the threshold to the highest auroc value in the eGCG_aucs vector if input 
  # threshold is too high 
  if(sum(eGCG_aucs >= auc_thres) == 0){
    warning(paste0('\n', 'The highest AUROC value between clusters for ', gene, 
    ' does not satisfy the input AUROC threshold, using the cluster with the highest AUROC value as the only eGCG+ cluster'))
    auc_thres <- max(eGCG_aucs)
  }
  
  # identify eGCG_negative and eGCG_positive clusters
  eGCG_neg<-names(eGCG_aucs[eGCG_aucs<auc_thres])
  eGCG_pos<-names(eGCG_aucs[eGCG_aucs>=auc_thres])
  
  # identify cells in eGCG_negative clusters and the least eGCG_positive cluster
  low_pos<-eGCG_pos[1]
  eGCG_neg_cells<-WhichCells(object,ident=eGCG_neg)
  eGCG_pos_cells<-WhichCells(object,ident=low_pos)
  
  # obtain count values
  obj_exp<-GetAssayData(object,slot = 'counts')
  eGCG_neg_cells_exps<-obj_exp[gene,eGCG_neg_cells]
  eGCG_pos_cells_exps<-obj_exp[gene,eGCG_pos_cells]
  
  # match the count value with group information and sort
  exp<-c(eGCG_neg_cells_exps,eGCG_pos_cells_exps)
  class<-c(rep(0,length(eGCG_neg_cells)),rep(1,length(eGCG_pos_cells)))
  names(exp)<-class
  exps<-sort(exp)
  
  # locate the best cutting-point
  quiet(curve<-roc(response = names(exps), predictor = exps, direction = '<'))
  yd_idx<-curve[["sensitivities"]]+curve[["specificities"]]-1
  roc_thres<-curve[["thresholds"]][which.max(yd_idx)]
  
  return(roc_thres)
}

#' Perform the contamination correction
#'
#' Decontaminate the count matrix of the input SeuratObject based on the input contaminative
#' genes using a Youden index-based method.
#'
#' @param object a clustered SeuratObject
#' @param cont_genes a contaminative geneset within the input SeuratObject
#' @param auc_thres the AUROC threshold to determine the boundary between eGCG_positive and eGCG_negative clusters (Default as 0.9, 90 percent)
#' @param min.cell the parameter used to filter the cell populations without sufficient number of cells. Cell populations that reaches the threshold could be used in downstream analysis.
#' 
#' @return the input clustered SeuratObject with a additional corrected assay of counts
#'
#' @import Seurat pROC ddpcr pbapply
#' @export
#'
ContaminationCorrection<-function(
    object,
    cont_genes,
    auc_thres = 0.9,
    min.cell = 50
){
  # calculate the threshold for every cont_gene (GCG)
  message('Calculating correction threshold...')
  object <- NormalizeData(object, normalization.method = "LogNormalize", 
                          scale.factor = 10000, verbose = F) # normalization
  num_cells<-table(Idents(object))
  qualified_cls<-names(num_cells[num_cells>=min.cell])
  thres_vals<-unlist(pblapply(cont_genes, function(x){
    ###
    quiet(eGCG_aucs<-Cal_AUCs(object,x,qualified_cls))
    thres<-Cal_thres(object,x,eGCG_aucs[[1]],auc_thres = auc_thres)
    return(thres)
  }))
  object@assays[[DefaultAssay(object)]]@data <- object@assays[[DefaultAssay(object)]]@counts # recover the status
  
  # fetch the matrix with cont_genes
  exp_matrix<-GetAssayData(object,slot='counts')
  decont_matrix_tmp<-exp_matrix[cont_genes,]
  if(length(cont_genes)==1){
    decont_matrix_tmp<-t(as.matrix(decont_matrix_tmp))
  }
  
  operation_matrix<-cbind(decont_matrix_tmp,thres=as.numeric(thres_vals))
  # correcting the counts
  message('Decontaminating...')
  corrected_mat_part<-pbapply(operation_matrix,1, function(x){
    thres<-x[length(x)]
    x<-x[1:(length(x)-1)]
    decont_x<-pmax(x-round(thres),0)
    return(decont_x)
  })
  corrected_mat_part<-t(as.matrix(corrected_mat_part))
  # cover the original counts with corrected counts
  exp_matrix[cont_genes,]<-corrected_mat_part
  
  Corrected_Assay<-CreateAssayObject(counts = Matrix(exp_matrix,sparse = T))
  object@assays[['Corrected']]<-Corrected_Assay
  # add the default key to avoid some version problems
  object@assays$Corrected@key = "corrected_"
  return(object)
}


#' Calculate contamination level using one GCG
#'
#' For one GCG, the contamination level is calculated by dividing its total expression in eGCG- cells by the total expression
#' of all genes in the eGCG- cells
#'
#' @param object a clustered SeuratObject
#' @param gene a gene within the input SeuratObject
#' @param eGCG_aucs a named vector of the AUROC values between the each cluster and the cluster with the lowest expression level
#' @param auc_thres the AUROC threshold to determine the boundary between eGCG_positive and eGCG_negative clusters (Default as 0.9, 90 percent)
#' @param slot the slot used for calculating contamination index
#' 
#' @return an index representing the contamination level calculated using the input gene
#'
#' @import Seurat 
#'
Cal_Cont_level <- function(object,gene,eGCG_aucs,auc_thres,slot){
  # setting the threshold to the highest auroc value in the eGCG_aucs vector if input 
  # threshold is too high 
  if(sum(eGCG_aucs >= auc_thres) == 0){
    auc_thres <- max(eGCG_aucs)
    # return(NA)
  }
  
  # identify eGCG_negative and eGCG_positive clusters
  eGCG_neg<-names(eGCG_aucs[eGCG_aucs<auc_thres])
  eGCG_pos<-names(eGCG_aucs[eGCG_aucs>=auc_thres])
  
  # identify cells in eGCG_negative clusters and the least eGCG_positive cluster
  # low_pos<-eGCG_pos[1]
  eGCG_neg_cells<-WhichCells(object,ident=eGCG_neg)
  eGCG_pos_cells<-WhichCells(object,ident=eGCG_pos)
  
  # obtain count values
  obj_exp<-GetAssayData(object,slot = slot)
  eGCG_neg_cells_exps<-obj_exp[gene,eGCG_neg_cells]
  eGCG_pos_cells_exps<-obj_exp[gene,eGCG_pos_cells]
  
  cont_index = log(sum(eGCG_neg_cells_exps)/sum(obj_exp[,eGCG_neg_cells]))+11
  
  return(cont_index)
}


#' Calculate contamination index for a dataset
#'
#' Calculating the contamination index for a dataset leveraging GCGs. For each GCG, the contamination level is calculated by dividing its total expression in eGCG- cells by 
#' the total expression of all genes in the eGCG- cells. The maximum of contamination index is returned. Higher contamination index represents higher contamination level.
#'
#' @param object a clustered SeuratObject
#' @param cont_genes a contaminative geneset within the input SeuratObject
#' @param min.cell the parameter used to filter the cell populations without sufficient number of cells. Cell populations that reaches the threshold could be used in downstream analysis.
#' @param top10_gcg a parameter controlling whether to use the top 10 GCGs for calculating contamination index
#' @param slot the slot used for calculating contamination index
#' 
#' @return the maximum value of the calculated contamination index
#'
#' @import Seurat pROC ddpcr pbapply
#' @export
#'
ContaminationQuantification <- function(object,cont_genes,
                         auc_thres = 0.9,
                         min.cell = 50,
                         top10_gcg = F,
                         slot = 'counts'){
  message('Calculating contamination level index...')
  object <- NormalizeData(object, normalization.method = "LogNormalize", 
                          scale.factor = 10000, verbose = F) # normalization
  num_cells<-table(Idents(object))
  if (top10_gcg ==T){
    if (length(cont_genes)>=10){
      cont_genes = cont_genes[1:10]
    }else
      cont_genes = cont_genes
  }
  qualified_cls<-names(num_cells[num_cells>=min.cell])
  index_vals<-unlist(pblapply(cont_genes, function(x){
    ###
    quiet(eGCG_aucs<-Cal_AUCs(object,x,qualified_cls))
    index<-Cal_Cont_level(object,x,eGCG_aucs[[1]],auc_thres = auc_thres,slot)
    return(index)
  }))
  if (max(index_vals)>4){
    message(paste0('The maximum contamination index is:',round(max(index_vals),2),', which indicates high contamination in this dataset. scCDC is highly recommended to be applied.'))
  }else{
    message(paste0('The maximum contamination index is:',round(max(index_vals),2),', which indicates low contamination in this dataset. Either DecontX or scCDC is recommended to be applied.'))
  }
  return(max(index_vals))
}