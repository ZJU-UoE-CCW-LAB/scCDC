% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Contamination_Correction.R
\name{ContaminationCorrection}
\alias{ContaminationCorrection}
\title{Perform the contamination correction}
\usage{
ContaminationCorrection(object, cont_genes, auc_thres = 0.9, min.cell = 50)
}
\arguments{
\item{object}{a clustered SeuratObject}

\item{cont_genes}{a contaminative geneset within the input SeuratObject}

\item{auc_thres}{the AUROC threshold to determine the boundary between eGCG_positive and eGCG_negative clusters (Default as 0.9, 90 percent)}

\item{min.cell}{the parameter used to filter the cell populations without sufficient number of cells. Cell populations that reaches the threshold could be used in downstream analysis.}
}
\value{
the input clustered SeuratObject with a additional corrected assay of counts
}
\description{
Decontaminate the count matrix of the input SeuratObject based on the input contaminative
genes using a Youden index-based method.
}
