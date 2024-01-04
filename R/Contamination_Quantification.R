#' Calculate contamination level using one GCG
#'
#' For one GCG, the contamination level is calculated by dividing its total expression in eGCG- cells by the total expression
#' of all genes in the eGCG- cells
#'
#' @param object a clustered SeuratObject
#' @param gene a gene within the input SeuratObject
#' @param eGCG_aucs a named vector of the AUROC values between the each cluster and the cluster with the lowest expression level
#' @param auc_thres the AUROC threshold to determine the boundary between eGCG_positive and eGCG_negative clusters (Default as 0.9, 90 percent)
#' @param slot the slot used for calculating contamination ratio
#' 
#' @return an ratio representing the contamination level calculated using the input gene
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
  
  cont_index = sum(eGCG_neg_cells_exps)/sum(obj_exp[,eGCG_neg_cells])
  
  return(cont_index)
}


#' Calculate contamination ratio for a dataset
#'
#' Calculating the contamination ratio for a dataset leveraging GCGs. For each GCG, the contamination level is calculated by dividing its total expression in eGCG- cells by 
#' the total expression of all genes in the eGCG- cells. The maximum of contamination ratio is returned. Higher contamination ratio represents higher contamination level.
#'
#' @param object a clustered SeuratObject
#' @param cont_genes a contaminative geneset within the input SeuratObject
#' @param min.cell the parameter used to filter the cell populations without sufficient number of cells. Cell populations that reaches the threshold could be used in downstream analysis.
#' @param top10_gcg a parameter controlling whether to use the top 10 GCGs for calculating contamination ratio
#' @param slot the slot used for calculating contamination ratio
#' 
#' @return the maximum value of the calculated contamination ratio
#'
#' @import Seurat pROC ddpcr pbapply
#' @export
#'
ContaminationQuantification <- function(object,cont_genes,
                                        auc_thres = 0.9,
                                        min.cell = 50,
                                        top10_gcg = F,
                                        slot = 'counts'){
  message('Calculating contamination ratio...')
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
  if (max(index_vals)>0.0003){
    message(paste0('The maximum contamination ratio is:',round(max(index_vals),2),', which indicates high contamination in this dataset. scCDC is highly recommended to be applied.'))
  }else{
    message(paste0('The maximum contamination ratio is:',round(max(index_vals),2),', which indicates low contamination in this dataset. Either DecontX or scCDC is recommended to be applied.'))
  }
  return(max(index_vals))
}
