## scCDC
library(scCDC)
set.seed(66)

data(mislet_before, package = "scCDC")
data(mislet_after, package = "scCDC")
test_sj <- CreateSeuratObject(mislet_before)
test_sj@active.ident <- mislet_annotation

test_that('ContaminationDetection with default settings works',{
  contamination <- ContaminationDetection(seuratobject = test_sj)
  expect_equal(contamination, mislet_cont)
  })

test_that('ContaminationDetection with default settings works',{
  corrected_mat <- ContaminationCorrection(object = test_sj, 
                                           cont_genes = rownames(mislet_cont))
  expect_equal(Matrix(corrected_mat, sparse = T), mislet_after)
})
