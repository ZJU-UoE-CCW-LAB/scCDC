# scCDC (single-cell Contamination Detection and Correction)
### Developed by Weijian Wang, Yihui Cen, Zezhen Lu

## Description
scCDC is a computational algorithm developed to detect global contamination-causing genes (GCGs) in single cell and single nuclei RNA-Seq datasets and perform further decontamination on the GCGs.

## Installation

`scCDC` can be installed from Github with the following code in `R`:

``` R
if(!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}

library(devtools)
install_github("ZJU-UoE-CCW-LAB/scCDC",ref = "sccdc-1.3")
```
This version of scCDC is developed based on Seurat V4. Seurat V4 can be installed with the following code in `R`:
``` R
if (!require("remotes", quietly = TRUE)){
    install.packages("remotes")
}

library(remotes)
remotes::install_github("satijalab/seurat", ref="release/4.3.0")
```
If you want to use Seurat V5, please check [scCDC version 1.4](https://github.com/ZJU-UoE-CCW-LAB/scCDC)
## Quick start
If you have a Seurat Object that contains clustering information, the typical scCDC workflow would be:


``` R
library(scCDC)
seuratobject = readRDS('/path/to/seuratobject')
GCGs = ContaminationDetection(seuratobject)
contamination_ratio = ContaminationQuantification(seuratobject,rownames(GCGs))
seuratobj_corrected = ContaminationCorrection(seuratobject,rownames(GCGs))
DefaultAssay(seuratobj_corrected) = "Corrected"
```

The decontaminated count matrix is stored in the 'Corrected' assay in the output Seurat Object, which can be directly used for downstram analysis. If you want to get the decontaminted count matrix, use the following code: 
```R
corrected_count_matrix = data.frame(seuratobj_corrected@assays[["Corrected"]]@counts)
```
If you want to start with count matrix, see [vignettes](https://github.com/ZJU-UoE-CCW-LAB/scCDC/tree/sccdc-1.3/blob/main/inst/doc/scCDC.html) for details.

## Usage

For detailed info on `scCDC` method and applications, please check out the package [vignettes](https://htmlpreview.github.io/?https://github.com/ZJU-UoE-CCW-LAB/scCDC/blob/main/inst/doc/scCDC.html), or with the following code in `R`: 

``` R
browseVignettes("scCDC")
```

## Contact

Any questions or suggestions on `scCDC` are welcomed! Please report it on [issues](https://github.com/ZJU-UoE-CCW-LAB/scCDC/issues), or contact Weijian Wang (<weijianwang@ucla.edu>), Yihui Cen (<yihuicen@g.ucla.edu>) or Zezhen Lu (<12307092@zju.edu.cn>).

## Reference
Wang, W., Cen, Y., Lu, Z., Xu, Y., Sun, T., Xiao, Y., Liu, W., Li, J. J., & Wang, C. (2024). scCDC: a computational method for gene-specific contamination detection and correction in single-cell and single-nucleus RNA-seq data. Genome biology, 25(1), 136. https://doi.org/10.1186/s13059-024-03284-w
