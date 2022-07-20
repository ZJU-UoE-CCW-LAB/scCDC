#' Calculate threshold using Otsu's method
#'
#' CalOtsu calculates proper threshold to seperate the distribution into 2 classes
#' using Otsu's method.
#' You can learn more about this function at:
#'
#' @param vec a numeric vector containing the whole distribution of counts
#'
#' @return a number indicating threshold calculated by Otsu's method
#'
#' @export
#'
#' @examples counts <- c(1, 1, 1, 2, 2, 2, 3, 3, 3);thres <- CalOtsu(counts)
# calculate the threshold using Otsu's method
CalOtsu <- function(vec){
  res_tmp <- c()
  data <- log(sort(vec) + 1)
  for(i in data){
    c1 <- data[which(data < i)]
    c2 <- data[which(data >= i)]
    v1 <- var(c1)
    v2 <- var(c2)
    p1 <- length(c1) / length(data)
    p2 <- length(c2) / length(data)
    tmp <- p1 * v1 + p2 * v2
    names(tmp) <- i
    res_tmp <- c(res_tmp, tmp)
  }
  res <- names(which.min(res_tmp))
  
  return(round(exp(as.numeric(res)) - 1))
}


#' Correction of contaminated counts based on CDF
#'
#' correspond_correction corrects the contaminated counts based on corresponding
#' contaminative counts distribution using a CDF-based method
#' You can learn more about this function at:
#'
#' @param target_cluster_exp numeric vector, the target counts distribution you want to decontaminate
#' @param cont_cluster_exp numeric vector, the corresponding contaminative counts distribution
#'
#' @return corrected vector of target count distribution
#'
#' @export
#'
#' @examples target <- c(5, 6, 7, 8);cont <- c(1, 2, 3, 4);corrected <- correspond_correction(target, cont)
# correction based on CDF
correspond_correction <- function(target_cluster_exp, cont_cluster_exp){
  target_sort <- sort(target_cluster_exp)
  cont_sort <- sort(cont_cluster_exp)
  corrected_list <- lapply(1:length(target_cluster_exp), function(index){
    
    # obtain the location(percentile) of each count from target cluster in target_sort 
    dense <- max(which(target_sort == target_cluster_exp[index])) / length(target_sort)
    # obtain the contamination count by the obtained location(percentile)
    correspond_cont <- cont_sort[ceiling(length(cont_cluster_exp) * dense)]
    
    return(max(target_cluster_exp[index] - correspond_cont, 0))
  })
  corrected <- unlist(corrected_list)
  names(corrected) <- names(target_cluster_exp)
  return(corrected)
}


#' Estimates the contaminative distribution
#'
#' Detecting_contaminated_cells find and expand the contaminated cluster to estimate the whole contamination distribution
#' You can learn more about this function at:
#'
#' @param x a named vector, values are normalized counts of a gene, names are the cell barcodes
#' @param cluster_info a named vector, names are the cell barcodes, values are the corresponding cluster of the cells,
#' usually are the active.ident in a Seurat object
#' @param was.thres a number indicating the threshold of wasserstein distance for
#' estimation of contaminative counts distribution
#'
#' @return a character vector containing the barcodes of cells that constitute 
#' the contamination distribution
#'
#' @import Seurat transport tibble dplyr pbapply
#' @export
#'
#' @examples cont_cells <- Detecting_contaminated_cells(x, Idents(object), 1)
# Determine the contaminative counts distribution
Detecting_contaminated_cells <- function(x, cluster_info, was.thres){
  total_clusters <- unique(cluster_info)
  # find the cluster with the lowest expression as the initial contamination
  exp_tmp <- data.frame(exp = x, cluster = cluster_info[names(x)])
  mid_mean_dat <- exp_tmp %>% group_by(cluster) %>% summarize(mid = median(exp), mean = mean(exp)) %>% arrange(mid, mean)
  cluster_sort <- as.vector(mid_mean_dat$'cluster')
  cont_cluster_out <- cluster_sort[1]
  
  if(cont_cluster_out == cluster_sort[length(cluster_sort)]){
    warning(paste0('No qualified cluster for a gene, no correction for this gene.'))
    cont_cells_expand <- NA
  }else{
    
    # get the cells in the lowest cluster 
    to_test_cells <- names(cluster_info)[cluster_info == cont_cluster_out]
    
    # calculate the wasserstein distance between the lowest cluster and other clusters
    was.list <- lapply(setdiff(total_clusters, cont_cluster_out), function(y){
      test_cells <- names(cluster_info[cluster_info == y])
      wasserstain  <-  wasserstein1d(as.numeric(x[test_cells]), as.numeric(x[to_test_cells]))
      return(wasserstain)
    })
    was.vals <- as.numeric(unlist(was.list))
    names(was.vals) <- setdiff(total_clusters, cont_cluster_out)
    # print(cont_cluster_out)
    # print(was.vals)
    # expand the distribution according to wasserstein distance
    cont_cluster_expand <- c(names(was.vals)[was.vals <= was.thres], cont_cluster_out)
    cont_cells_expand <- names(cluster_info)[cluster_info %in% cont_cluster_expand]
    names(cont_cells_expand) <- as.character(cluster_info[cont_cells_expand])
  }
  # return the the cell barcodes
  return(cont_cells_expand)
}



#' Corrects the contamination in Seurat
#'
#' DecontFunc estimates the contamination of the given genes and corrects them using
#' Otsu's method or CDF-based method
#' You can learn more about this function at:
#'
#' @param object a Seurat object that has been clustered
#' @param cont_genes a character vector, containing the genes to be corrected, 
#' usually generated by ContaminationDetection
#' @param was.thres a number indicating the threshold of wasserstein distance for
#' estimation of containative counts distribution
#'
#' @return a decontaminated count matrix 
#'
#' @import Seurat tibble pbapply
#' @export
#'
#' @examples cont_genes <- c('Ins1', 'Ins2');DecontFunc(object, cont_genes, 1)
DecontFunc <- function(
  object,
  cont_genes,
  was.thres
){
  # PREPARATION
  # extract count matrix
  exp_matrix <- as.matrix(GetAssayData(object, slot = 'counts'))
  
  # generate new normalized count matrix
  cell_factors <- 10000 / colSums(exp_matrix)
  exp_matrix_list <- lapply(1:ncol(exp_matrix), function(x){
    return(exp_matrix[, x] * cell_factors[x])
  })
  exp_matrix_norm <- do.call(cbind, exp_matrix_list)
  rownames(exp_matrix_norm) <- rownames(object)
  colnames(exp_matrix_norm) <- colnames(object)
  exp_matrix_norm <- log(exp_matrix_norm + 1)
  
  # get cluster information
  cluster_info <- Idents(object)
  total_clusters <- unique(cluster_info)
  
  # extract count and normalized count matrix for only contaminated genes entered
  if(length(cont_genes) == 1){
    exp_matrix_cont <- t(as.matrix(exp_matrix[cont_genes, ]))
    exp_matrix_cont_norm <- t(as.matrix(exp_matrix_norm[cont_genes, ]))
  }else{
    exp_matrix_cont <- exp_matrix[cont_genes, ]
    exp_matrix_cont_norm <- exp_matrix_norm[cont_genes, ]
  }
  
  #CORRECTION
  message('Decontaminating...')
  
  # select fully contaminted cells for contamination estimation (normalized count matrix used)
  cont_cells_tmp <- pbapply(exp_matrix_cont_norm, 1, function(x){
    cont_cells <- Detecting_contaminated_cells(x, cluster_info, was.thres)
    return(cont_cells)
  })
  if(is.null(dim(cont_cells_tmp))){
    cont_cells_mat <- as.matrix(cont_cells_tmp)
  }else{
    cont_cells_mat <- t(cont_cells_tmp)
  }
  # generate the corrected matrix for contaminated genes entered (count matrix used)
  gene_index <- 1:length(cont_genes)
  exp_matrix_cont <- cbind(exp_matrix_cont, gene_index)
  corrected_exp_mat <- apply(exp_matrix_cont, 1, function(x){
    gene_tmp <- x['gene_index']
    x <- x[1:(length(x) - 1)]
    cell_arrange <- names(x)
    cont_cells <- unlist(cont_cells_mat[gene_tmp, ])
    
    
    ### judge if otsu should be used
    cont_cls <- unique(gsub('.*\\.', '', names(cont_cells)))
    if(length(setdiff(total_clusters, cont_cls)) == 1){ # when only one cluster is special
      # print(c(gene_tmp, 'Otsu'))
      dif_cluster <- setdiff(total_clusters, cont_cls)
      dif_cluster_cells <- WhichCells(object, idents = dif_cluster)
      
      # bootstrap to make the cell number in each cluster equal
      boot_num <- max(length(cont_cells), length(dif_cluster_cells))
      cont_boot <- sample(x[cont_cells], boot_num, replace = T)
      dif_cluster_boot <- sample(x[dif_cluster_cells], boot_num, replace = T)
      
      # combined the distribution to apply Otsu's method
      boot_combined <- c(cont_boot, dif_cluster_boot)
      #plot(hist(log(boot_combined+1), breaks=1000))
      Otsu_thres <- CalOtsu(boot_combined)
      # print(Otsu_thres)
      corrected_tmp <- x - Otsu_thres
      corrected_tmp[which(corrected_tmp < 0)] <- 0
      corrected_vec <- corrected_tmp
      return(corrected_vec)
    }
    ###
    
    
    
    # otherwise use the CDF-based correction
    correction_list <- lapply(total_clusters, function(y){
      decont_cells <- WhichCells(object, idents = y)
      target_cluster_exp <- x[decont_cells]
      cont_cluster_exp <- x[cont_cells]
      #print(c(gene_tmp, 'corres'))
      corrected <- correspond_correction(target_cluster_exp, cont_cluster_exp)
      return(corrected)
    })
    corrected_vec <- unlist(correction_list)[cell_arrange]
    return(corrected_vec)
  })
  corrected_exp_mat <- t(corrected_exp_mat)
  
  # replace the observed counts by corrected counts of contaminated genes
  exp_matrix[cont_genes, ] <- corrected_exp_mat
  
  message('Done!')
  
  return(exp_matrix)
}  

#' Interface to use
#'
#' ContaminationCorrection is the interface function to correct the contamination
#' in single-cell sequencing data
#' You can learn more about this function at:
#'
#' @param object a Seurat object that has been clustered
#' @param cont_genes a character vector, containing the genes to be corrected, 
#' usually generated by ContaminationDetection
#' @param was.thres a number indicating the threshold of wasserstein distance for
#' estimation of containative counts distribution, default is 1
#'
#' @return a decontaminated count matrix
#'
#' @export
#'
#' @examples cont_genes <- rownames(ContaminationDetection(object));DecontFunc(object, cont_genes, 1)
#INTERFACE
ContaminationCorrection <- function(
  object,
  cont_genes,
  was.thres = 1
){
  decont_matrix <- DecontFunc(object, cont_genes, was.thres)
  return(decont_matrix)
}