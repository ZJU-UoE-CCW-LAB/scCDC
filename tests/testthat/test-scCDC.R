## scCDC
library(scCDC)
library(Seurat)
set.seed(66)

data(mislet_before, package = "scCDC")
data(mislet_after, package = "scCDC")
test_sj <- CreateSeuratObject(mislet_before)
test_sj@active.ident <- mislet_annotation

test_that('ContaminationDetection with default settings works',{
  contamination <- ContaminationDetection(seuratobject = test_sj, 
                                          restriction_factor = 0.8)
  expect_equal(contamination, mislet_cont)
  })

test_that('ContaminationCorrection with default settings works',{
  test_sj_corrected <- ContaminationCorrection(object = test_sj, 
                                           cont_genes = rownames(mislet_cont),
                                           auc_thres = 0.9)
  DefaultAssay(test_sj_corrected) <- "Corrected"
  test_sj_corrected_mat <- GetAssayData(test_sj_corrected, assay = "Corrected", 
                                        slot = "counts")
  expect_equal(test_sj_corrected_mat, mislet_after)
})
