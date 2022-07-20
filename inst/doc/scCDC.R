## ----echo=FALSE, include=FALSE------------------------------------------------
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, dev = "png",
                      message = FALSE, error = FALSE, warning = FALSE)

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(scCDC))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(pdftools))
suppressPackageStartupMessages(library(ddpcr))

## ----load-data----------------------------------------------------------------
data(mislet_before, package = "scCDC")
mislet_seuratobj<-CreateSeuratObject(mislet_before)
mislet_seuratobj@active.ident<-mislet_annotation

## ----ContaminationDetection.1-------------------------------------------------
contamination <- ContaminationDetection(mislet_seuratobj,sample_name = "mislet",
                                      out_path.plot = "./", out_path.table = "./")
rownames(contamination)

## ----ContaminationDetection.2,fig.width=1,fig.height=1------------------------
quiet(pdf_convert("./mislet_SE-plot.pdf",format='jpg'))
knitr::include_graphics("mislet_SE-plot_1.jpg")
knitr::include_graphics("mislet_SE-plot_2.jpg")
knitr::include_graphics("mislet_SE-plot_3.jpg")
knitr::include_graphics("mislet_SE-plot_4.jpg")

## ----ContaminationDetection.3-------------------------------------------------
cont_table <- read.csv("mislet_degree_of_contamination.csv")
cont_table_long <- tidyr::gather(cont_table, group, entropy_deviation, 
                                 Alpha, Beta, Delta, Endothelial, mean_distance)
cont_table_long_selected <- cbind(cont_table_long, 
                                  type = rep("uncontaminated", length(cont_table_long[,1])))
cont_table_long_selected[which(cont_table_long_selected$X %in% 
                                 rownames(contamination)),4] <- "contaminated"
ggplot(cont_table_long_selected) + geom_boxplot(aes(type,entropy_deviation),outlier.shape = NA) +
  facet_wrap(~group)


## ----ContaminationCorrection.1------------------------------------------------
mislet_corrected <- ContaminationCorrection(mislet_seuratobj, 
                                                      cont_genes = rownames(contamination))

## ----ContaminationCorrection.2,message=FALSE,warning=FALSE--------------------
# before the decontamination
mislet_seuratobj <- NormalizeData(mislet_seuratobj,
                                  normalization.method = "LogNormalize",scale.factor = 10000)
mislet_seuratobj <- FindVariableFeatures(mislet_seuratobj, 
                                              selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mislet_seuratobj)
mislet_seuratobj <- ScaleData(mislet_seuratobj, features = all.genes)
mislet_seuratobj <- RunPCA(mislet_seuratobj, 
                                features = VariableFeatures(object = mislet_seuratobj))

mislet_seuratobj <- RunUMAP(mislet_seuratobj,dims=1:20)
# after the decontamination
mislet_seuratobj_corrected <- CreateSeuratObject(mislet_corrected)
mislet_seuratobj_corrected <- NormalizeData(mislet_seuratobj_corrected,
                                  normalization.method = "LogNormalize",scale.factor = 10000)
mislet_seuratobj_corrected <- FindVariableFeatures(mislet_seuratobj_corrected, 
                                              selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mislet_seuratobj_corrected)
mislet_seuratobj_corrected <- ScaleData(mislet_seuratobj_corrected, features = all.genes)
mislet_seuratobj_corrected <- RunPCA(mislet_seuratobj_corrected, 
                                features = VariableFeatures(object = mislet_seuratobj_corrected))

mislet_seuratobj_corrected <- RunUMAP(mislet_seuratobj_corrected,dims=1:20)

## ----ContaminationCorrection.3------------------------------------------------
FeaturePlot(mislet_seuratobj, features = c("Ins1", "Gcg"))
FeaturePlot(mislet_seuratobj_corrected, features = c("Ins1", "Gcg"))

## -----------------------------------------------------------------------------
sessionInfo()
