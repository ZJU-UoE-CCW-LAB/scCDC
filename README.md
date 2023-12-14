# scCDC (single-cell Contamination Detection and Correction)
### Developed by Weijian Wang, Yihui Cen, Zezhen Lu

## Description
scCDC is a computational algorithm developed to detect global contamination causing genes (GCGs) in single cell and single nuclei RNA-Seq datasets and perform further decontamination on the GCGs.

## Installation

`scCDC` can be installed from Github with the following code in `R`:

``` r
if(!require("devtools")){
  install.packages("devtools")
}

library(devtools)
install_github("ZJU-UoE-CCW-LAB/scCDC")
```
## Quick start
If you have a Seurat Object that contains clustering information, the typical scCDC workflow would be:

```{r Quick_start,eval=FALSE}
# load data
seuratobject = readRDS('/path/to/seuratobject')
# detect global contamination causing genes(GCGs)
GCGs <- ContaminationDetection(seuratobject)
# remove the contamination
seuratobj_corrected <- ContaminationCorrection(seuratobject,rownames(GCGs))
DefaultAssay(seuratobj_corrected) <- "Corrected"
```

The decontaminated count matrix is stored in the 'Corrected' assay in the output Seurat Object, which can be directly used for downstram analysis. If you want to get the decontaminted count matrix, use the following code: 
```{r eval=FALSE}
corrected_count_matrix = data.frame(seuratobj_corrected@assays[["Corrected"]]@counts)
```

## Usage

For detailed info on `scCDC` method and applications, please check out the package [vignettes](https://htmlpreview.github.io/?https://github.com/ZJU-UoE-CCW-LAB/scCDC/blob/main/inst/doc/scCDC.html), or with the following code in `R`: 

``` r
browseVignettes("scCDC")
```

## Contact

Any questions or suggestions on `scCDC` are welcomed! Please report it on [issues](https://github.com/ZJU-UoE-CCW-LAB/scCDC/issues), or contact Weijian Wang (<weijian.19@intl.zju.edu.cn>), Yihui Cen (<yihui.19@intl.zju.edu.cn>) or Zezhen Lu (<zezhen.19@intl.zju.edu.cn>).

## Reference
