% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Contamination_Detection.R
\name{CalExpressionPercent.Seurat}
\alias{CalExpressionPercent.Seurat}
\title{Calculate expression percent}
\usage{
CalExpressionPercent.Seurat(object, gene_set)
}
\arguments{
\item{object}{a Seurat object that has been clustered.}

\item{gene_set}{a vector containing genes which you want to calculate the expression percent.}
}
\value{
a matrix containing expression percent results with gene names as row names, cluster names as column names.
}
\description{
CalExpressionPercent.Seurat calculates the expression percent of genes within each cluster 
based on the 'counts' slot of the given seurat object.
}
