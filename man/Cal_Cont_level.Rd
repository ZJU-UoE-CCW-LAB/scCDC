% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Contamination_Quantification.R
\name{Cal_Cont_level}
\alias{Cal_Cont_level}
\title{Calculate contamination level using one GCG}
\usage{
Cal_Cont_level(object, gene, eGCG_aucs, auc_thres, slot)
}
\arguments{
\item{object}{a clustered SeuratObject}

\item{gene}{a gene within the input SeuratObject}

\item{eGCG_aucs}{a named vector of the AUROC values between the each cluster and the cluster with the lowest expression level}

\item{auc_thres}{the AUROC threshold to determine the boundary between eGCG_positive and eGCG_negative clusters (Default as 0.9, 90 percent)}

\item{slot}{the slot used for calculating contamination ratio}
}
\value{
an ratio representing the contamination level calculated using the input gene
}
\description{
For one GCG, the contamination level is calculated by dividing its total expression in eGCG- cells by the total expression
of all genes in the eGCG- cells
}
