## ----echo=FALSE, include=FALSE------------------------------------------------
library(knitr)
knitr::opts_chunk$set(tidy = FALSE, dev = "png",
                      message = FALSE, warning = FALSE)

## ----setup,echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages(library(scCDC))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ddpcr))

## ----Quick_start,eval=FALSE---------------------------------------------------
#  library(devtools)
#  install_github("ZJU-UoE-CCW-LAB/scCDC")
#  library(scCDC)
#  # load data
#  seuratobject <- readRDS('/path/to/seuratobject')
#  # detect global contamination causing genes(GCGs)
#  GCGs <- ContaminationDetection(seuratobject)
#  # remove the contamination
#  seuratobj_corrected <- ContaminationCorrection(seuratobject,rownames(GCGs))
#  DefaultAssay(seuratobj_corrected) <- "Corrected"

## ----eval=FALSE---------------------------------------------------------------
#  corrected_count_matrix <- data.frame(seuratobj_corrected@assays[["Corrected"]]@counts)

## ----eval=FALSE---------------------------------------------------------------
#  if(!require("devtools")){
#    install.packages("devtools")
#  }
#  library(devtools)
#  install_github("ZJU-UoE-CCW-LAB/scCDC")

## -----------------------------------------------------------------------------
library(scCDC)

## ----load-data----------------------------------------------------------------
data(mislet_before, package = "scCDC")

## ----apply-Seurat,message=FALSE,warning=FALSE---------------------------------
mislet_seuratobj<-CreateSeuratObject(mislet_before)
mislet_seuratobj <- NormalizeData(mislet_seuratobj,
                                  normalization.method = "LogNormalize",scale.factor = 10000)
mislet_seuratobj <- FindVariableFeatures(mislet_seuratobj, 
                                              selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mislet_seuratobj)
mislet_seuratobj <- ScaleData(mislet_seuratobj, features = all.genes)
mislet_seuratobj <- RunPCA(mislet_seuratobj, 
                                features = VariableFeatures(object = mislet_seuratobj))
# We use pre-defined clustering annotation here
mislet_seuratobj@active.ident<-mislet_annotation
mislet_seuratobj <- RunUMAP(mislet_seuratobj,dims=1:20)

## ----ContaminationDetection.0-------------------------------------------------
GCGs <- ContaminationDetection(mislet_seuratobj)
rownames(GCGs)

## ----ContaminationDetection.1, eval=FALSE-------------------------------------
#  GCGs <- ContaminationDetection(mislet_seuratobj,restriction_factor = 0.5,
#                                          sample_name = "mislet",out_path.plot = "./",
#                                          out_path.table = "./")

## -----------------------------------------------------------------------------
mislet_cont_ratio <- ContaminationQuantification(mislet_seuratobj,rownames(GCGs))
mislet_cont_ratio

## ----ContaminationCorrection.1------------------------------------------------
mislet_seuratobj_corrected <- ContaminationCorrection(mislet_seuratobj, rownames(GCGs))

## ----eval=FALSE---------------------------------------------------------------
#  corrected_count_matrix = data.frame(mislet_seuratobj_corrected@assays[["Corrected"]]@counts)

## ----ContaminationCorrection.2------------------------------------------------
DefaultAssay(mislet_seuratobj_corrected) <- "Corrected"

## ----ContaminationCorrection.3,result='hide',fig.show='hide'------------------
mislet_seuratobj_corrected <- NormalizeData(mislet_seuratobj_corrected,
                                  normalization.method = "LogNormalize",scale.factor = 10000)
mislet_seuratobj_corrected <- FindVariableFeatures(mislet_seuratobj_corrected,
                                                   selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mislet_seuratobj_corrected)
mislet_seuratobj_corrected <- ScaleData(mislet_seuratobj_corrected,
                                        features = all.genes)
mislet_seuratobj_corrected <- RunPCA(mislet_seuratobj_corrected, 
                                features = VariableFeatures(object = mislet_seuratobj_corrected))
ElbowPlot(mislet_seuratobj_corrected)
mislet_seuratobj_corrected <- FindNeighbors(mislet_seuratobj_corrected, dims = 1:15)
mislet_seuratobj_corrected <- FindClusters(mislet_seuratobj_corrected, resolution = 0.2,verbose = F)
mislet_seuratobj_corrected <- RunUMAP(mislet_seuratobj_corrected,dims=1:15)

## ----ContaminationCorrection.4------------------------------------------------
FeaturePlot(mislet_seuratobj, features = c("Ins1", "Gcg"))
FeaturePlot(mislet_seuratobj_corrected, features = c("Ins1", "Gcg"))

## -----------------------------------------------------------------------------
sessionInfo()

