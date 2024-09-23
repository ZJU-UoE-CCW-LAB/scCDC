#' Calculate expression percent
#'
#' CalExpressionPercent.Seurat calculates the expression percent of genes within each cluster 
#' based on the 'counts' layer of the given seurat object.
#'
#' @param object a Seurat object that has been clustered.
#' @param gene_set a vector containing genes which you want to calculate the expression percent.
#' 
#' @return a matrix containing expression percent results with gene names as row names, cluster names as column names.
#'
#' @import Seurat
#'
CalExpressionPercent.Seurat <- function(
  object,
  gene_set
){
  # fetch the cluster information
  levels <- levels(object@active.ident)
  # calculate the expression percentage of each cluster for each gene
  output.list <- lapply(levels, function(x){
    tmp <- subset(object, idents = x)
    gene.data.dat <- as.data.frame(GetAssayData(tmp, layer = 'counts'))[gene_set, ]
    gene.exp.percent <- apply(gene.data.dat, 1, function(x){
      gene.exp.percent <- 1-(length(which(x == 0))/length(x))
      return(gene.exp.percent)
    })
    return(gene.exp.percent)
  })
  # reformat the result
  output <- do.call(cbind, output.list)
  colnames(output) <- levels
  return(output)
}

#' Calculate average expression
#'
#' CalAverageExpression.Seurat calculates average expression of genes within each cluster 
#' based on the data/counts layer of the given seurat object.
#' @param object a Seurat object that has been clustered.
#' @param gene_set a vector containing genes which you want to calculate average expression.
#' @param type a parameter setting for whether to use 'data' layer or 'counts' layer in the Seurat object to calculate average expression. "default" for calculation using 'counts' layer.
#'
#' @return a matrix containing average expression results with gene names as row names, cluster names as column names.
#'
#' @import Seurat
#' @import pbapply
#'
CalAverageExpression.Seurat <- function(
  object,
  gene_set,
  type = 'counts'
){
  # fetch the cluster information
  levels <- levels(object@active.ident)
  # calculate the mean expression level of each cluster for each gene
  output.list <- pblapply(levels, function(x) {
    tmp <- subset(object, idents = x)
    gene.data.dat <- as.data.frame(GetAssayData(tmp, layer = type))[gene_set, ]
    gene.set.data <- apply(gene.data.dat, 1, mean)
    return(gene.set.data)
  })
  # reformat the result
  output <- do.call(cbind, output.list)
  colnames(output) <- levels
  return(output)
}

#' Calculate the entropy 
#'
#' This method calculates the Shannon's entropy of the observed counts of genes within each cluster in the Seurat object.
#' CalEnt.Seurat calculates entropy of genes based on the 'counts' layer in the Seurat object.
#' 
#' @param object a Seurat object that has been clustered.
#' @param gene_set a vector containing genes which you want to calculate entropy for.
#' 
#' @return a matrix containing entropy results with gene names as row names, cluster names as column names.
#' 
#' @import Rcpp
#' @import Seurat
#' @import pbapply
#'
CalEnt.Seurat <- function(
  object,
  gene_set
){
  # fetch the cluster information
  levels <- levels(object@active.ident)
  # calculate the entropy of each cluster for each gene
  output.list <- pblapply(levels, function(x) {
    tmp <- subset(object, idents = x)
    gene.data.dat <- as.matrix(GetAssayData(tmp, layer = 'counts'))[gene_set, ]
    if(is.matrix(gene.data.dat)){
      gene.set.ent <- MatrixToEntropy(gene.data.dat)
    }
    else if(is.vector(gene.data.dat)){
      gene.set.ent <- VectorToEntropy(gene.data.dat)
    }else{
      stop("gene.data.dat is not a vector or matrix")
    }
    return(gene.set.ent)
  })
  # reformat the result
  output <- do.call(cbind, output.list)
  colnames(output) <- levels
  rownames(output) <- gene_set
  return(output)
}

#' Fit a curve representing the relationship between the entropy and the mean expression level.
#' 
#' generate_curve fits the relationship between entropy and mean expression level through bootstrapping
#' generate_curve estimates the expected entropy based on the mean expression level within each cluster in the Seurat object
#' 
#' @param .x a tibble object containing 3 columns: gene name, mean expression and entropy.
#' @param filter the parameter used to filter the points that are deviated from the curve.
#' @param select_factor the parameter controls the selected number of genes used to do bootstrap.
#' 
#' @return a dataframe containing 7 columns: mean expression, the actual entropy, entropy divergence, the fitted entropy, p value, the adjusted p value, gene name.
#' @import dplyr
#'
generate_curve <- function(.x, filter = 0.01, select_factor = 0.8){
  number <- round(nrow(.x)*select_factor)
  max_number <- max(.x$mean.expr)
  discrete_x <- seq(0, max_number+0.1, 0.001)
  ## every round use 80% genes to fit the curve, and then average these curves
  result <- lapply(1:10, function(i){
    tmp <- .x[sample(1:nrow(.x), number, replace = F), ]
    fit <- smooth.spline(tmp$mean.expr, tmp$entropy, spar = 1)
    prd <- predict(fit, tmp$mean.expr)
    tmp %>%
      dplyr::mutate(fit = prd[[2]]) %>%
      dplyr::mutate(distance = fit - entropy) %>%
      dplyr::filter(is.finite(distance)) %>%
      dplyr::mutate(adjusted_p = 1-pnorm(.$distance, mean = mean(.$distance), sd = sd(.$distance))) %>%
      dplyr::filter(adjusted_p > filter) -> tmp
    fit <- smooth.spline(tmp$mean.expr, tmp$entropy, spar = 1)
    prd <- predict(fit, tmp$mean.expr)
    tmp %>%
      dplyr::mutate(fit = prd[[2]]) %>%
      dplyr::mutate(distance = fit - entropy) %>%
      dplyr::filter(is.finite(distance)) %>%
      dplyr::mutate(adjusted_p = 1-pnorm(.$distance, mean = mean(.$distance), sd = sd(.$distance))) %>%
      dplyr::filter(adjusted_p > filter) -> tmp
    fit <- smooth.spline(tmp$mean.expr, tmp$entropy, spar = 1)
    prd <- predict(fit, discrete_x)
    return(prd[[2]])
  })
  result <- data.frame(do.call(cbind, result))
  result$fit <- rowMeans(result, na.rm = T)
  rownames(result) <- as.numeric(discrete_x)
  result <- result %>%
    dplyr::select(fit)
  
  mean.expr = as.numeric(round(.x$mean.expr, 3))
  entropy =.x$entropy
  tmp_fit <- result[as.character(mean.expr), 'fit']
  tmp <- data.frame(mean.expr = mean.expr, entropy = entropy, distance = tmp_fit-entropy, fit = tmp_fit)
  rownames(tmp) <- .x$Gene
  
  tmp <- tmp %>% dplyr::mutate(p.value = 1-pnorm(tmp$distance, mean = mean(tmp$distance), sd = sd(tmp$distance)))
  p.adj <- p.adjust(tmp$p.value, method = "fdr")
  tmp <- tmp %>% dplyr::mutate(p.adj = p.adj) %>% dplyr::arrange(desc(distance))
  tmp$Gene <- rownames(tmp)
  return(tmp)
}


#' Plot the relationship between the entropy and the expression
#' 
#' generate_plot plots the fitted curve with points. Each point represents a gene in a cluster, the x axis represents the mean expression level
#' and the y axis represents the entropy.
#' 
#' @param .x a dataframe which is the output of the generate_curve function.
#' @param name the name printed on the figure.
#' @param genes if specified a vector of genes, only these genes would be labeled on the plot. If not specified, the significantly deviated genes in the dataset would be labeled.
#' @param point_size the parameter controls the point size in the plot.
#' @param cutoff the threshold controls the justification of significantly deviated genes.
#'
#' @return A ggplot object
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#'
generate_plot <- function(.x, name, genes = NULL, point_size = 1.8, cutoff = 0.05){
  if(is.null(genes)){
    .x <- .x %>% dplyr::mutate(sig = ifelse(p.adj <= cutoff, 1, 0))
    con <- .x %>%
      filter(sig == 1)
    uncon <- .x %>%
      filter(sig == 0)
    p =
      ggplot()+
      geom_point(data = con, aes(mean.expr, entropy), alpha = 0.8, color = "#d75427", size = point_size) +
      geom_point(data = uncon, aes(mean.expr, entropy), alpha = 0.8, color = "#2e409a", size = point_size)+
      geom_text(data = con, aes(mean.expr, entropy, label = Gene), nudge_x = 0.3, nudge_y = 0, check_overlap = T)+
      geom_line(data = .x, aes(mean.expr, fit), lwd = 0.7) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 0),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.grid = element_blank()
      ) +
      labs(
        title = name,
        x = "ln(mean expression)",
        y = "Shannon's entropy of count"
      )+theme(plot.title = element_text(hjust = 0.5))
    return(p)
  }else{
    .x <- .x %>% dplyr::mutate(sig = ifelse(p.adj <= cutoff, 1, 0))
    con_data <- .x[genes, ]
    con <- .x %>%
      filter(sig == 1)
    uncon <- .x %>%
      filter(sig == 0)
    p =
      ggplot()+
      geom_point(data = uncon, aes(mean.expr, entropy), alpha = 0.8, color = "#2e409a", size = point_size) +
      geom_point(data = con, aes(mean.expr, entropy), alpha = 0.8, color = "#d75427", size = point_size)+
      geom_text_repel(data = con_data, aes(mean.expr, entropy, label = Gene), max.overlaps = 100)+
      geom_line(data = .x, aes(mean.expr, fit), lwd = 0.7) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 0),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.grid = element_blank()
      ) +
      labs(
        title = name,
        x = "ln(mean expression)",
        y = "Shannon's entropy of count"
      )+theme(plot.title = element_text(hjust = 0.5))
    return(p)
  }
}

#' Detect the comtamination causing genes 
#'
#' ContaminationDetection detects the contamination causing genes in the Seurat object based on the 
#' relationship between the entropy and mean expression level of genes in a cell cluster.
#'
#' @param seuratobject a Seurat object that has been clustered.
#' @param restriction_factor the parameter controls the degree of conservation when justifying the contamination causing genes. Default setting is 0.8, representing that each potential contamination-causing gene should be recognized in 80 percent of the clusters with the algorithm to be identified in the output result.
#' @param sample_name the name of the output contamination degree dataframe.
#' @param min.cell the parameter used to filter the cell populations without sufficient number of cells. Cell populations that reaches the threshold could be used in downstream analysis.
#' @param percent.cutoff the parameter used to filter the candidate contamination causing genes without sufficient expression percentage in each cluster
#' @param out_path.plot if specified a path, the plot of the relationship between the entropy and expression would be output into the path.
#' @param out_path.table if specified a path, a contamination degree dataframe would be output into the path.
#' @param only.cont_genes a logical parameter controls whether only to label the contamination causing genes on the entropy-expression curve
#' 
#' @return a dataframe containing entropy divergence results with gene names as row names, cluster names as column names, the last column is the mean entropy divergence of one gene across the clusters.
#'
#' @import ddpcr
#' @import pbapply
#' @import dplyr
#' @import Seurat
#' @export
#'
ContaminationDetection <- function(seuratobject, restriction_factor = 0.5, sample_name = "default",
                                         min.cell = 100, percent.cutoff = 0.2, out_path.plot = NULL, out_path.table = NULL, only.cont_genes = F){
  # fetch the the cell number of each cluster
  cluster <- table(seuratobject@active.ident)
  nums <- as.numeric(cluster)
  ##the number of cells in a cluster should reach a threshold
  total_cluster <- names(cluster)[nums >= min.cell]
  seuratobject <- subset(seuratobject, idents = total_cluster)
  
  # calculate entropy
  message("Caculating entropy...")
  entropy_result <- CalEnt.Seurat(seuratobject, rownames(seuratobject))
  # calculate average expression level
  message("Caculating expression level...")
  ave_result <- CalAverageExpression.Seurat(seuratobject, rownames(seuratobject))
  ave_result <- log(ave_result+1)
  total_cluster <- colnames(ave_result)
  # calculate the relation betwee entropy and expression level
  message("Calculating entropy-expression relation...")
  all <- pblapply(1:length(total_cluster), function(i){
    cluster <- total_cluster[i]
    gene_ent_ave <- tibble(Gene = rownames(entropy_result), mean.expr = ave_result[, which(colnames(ave_result) == cluster)], entropy = entropy_result[, which(colnames(entropy_result) == cluster)])
    gene_ent_ave <- generate_curve(gene_ent_ave)
    return(gene_ent_ave)
  })
  
  # select contaminated genes (GCG)
  genes <- c()
  for (i in 1:length(all)){
    tmp = all[[i]]
    genetmp <- tmp %>% filter(p.adj <= 0.05) %>%dplyr::select(Gene)
    genes <- c(genes, as.vector(as.matrix(genetmp)))
  }
  count <- table(genes)
  gene_num <- as.numeric(count)
  contaminated_genes <- names(count)[gene_num >= round(ncol(ave_result) * restriction_factor)]
  if (length(contaminated_genes)==0){
    stop('No contaminated genes found')
  }
  exp_percent <- CalExpressionPercent.Seurat(seuratobject, contaminated_genes)
  
  # filter contaminated genes (GCG)
  filtered_contaminated_genes <- c()
  for(i in 1:length(contaminated_genes)){
    choose <- T
    for (j in 1:ncol(exp_percent)){
      if(exp_percent[i, j]<percent.cutoff){
        choose <- F
      }
    }
    if (choose == T){
      filtered_contaminated_genes <- c(filtered_contaminated_genes, contaminated_genes[i])
    }
  }
  if (length(filtered_contaminated_genes)==0){
    stop('No contaminated genes found')
  }
  # extract cont_degree
  message("Extracting contamination degree...")
  genes <- rownames(seuratobject)
  distance_list <- pblapply(1:length(all), function(x){
    tmp <- all[[x]]
    distance <- tmp[genes, 'distance']
    return(distance)
  })
  distance_result <- as.data.frame(do.call(cbind, distance_list))
  colnames(distance_result) <- total_cluster
  rownames(distance_result) <- genes
  distance_result <- cbind(distance_result, mean_distance = rowMeans(distance_result))
  contamination_result <- distance_result %>%
    arrange(desc(as.numeric(mean_distance)))
  
  # output degree of contamination
  if (!is.null(out_path.table)){
    write.csv(contamination_result, paste0(out_path.table, sample_name, "_degree_of_contamination.csv"))
  }
    
  # output Entropy-Expression relation
  if (!is.null(out_path.plot)){
    if(only.cont_genes == F){
      pdf(paste0(out_path.plot, sample_name, '_SE-plot.pdf'))
      plot_list <- c()
      for (i in 1:length(all)){
        plot_list[i] <- list(generate_plot(all[[i]], total_cluster[i]))
      }
      for (i in 1:length(all)){
        quiet(print(plot_list[i]))
      }
      dev.off()
    }else{
      pdf(paste0(out_path.plot, sample_name, '_SE-plot.pdf'))
      plot_list <- c()
      for (i in 1:length(all)){
        plot_list[i] <- list(generate_plot(all[[i]], total_cluster[i], genes = filtered_contaminated_genes))
      }
      for (i in 1:length(all)){
        quiet(print(plot_list[i]))
      }
      dev.off()
    }
  }
   
  message(paste0("Complete detection. ", length(filtered_contaminated_genes), " contaminated genes found"))
  
  return(contamination_result[which(rownames(contamination_result) %in% filtered_contaminated_genes), ])
}