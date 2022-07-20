# scCDC
### Developed by Zezhen Lu, Weijian Wang, Yihui Cen

## Description
scCDC is a computational algorithm developed to detect contamination causing genes (GCGs) in single cell and single nuclei RNA-Seq datasets and perform further decontamination on the GCGs.

## Installation

`scCDC` can be installed from Github with the following code in `R`:

``` r
install.packages("devtools")
library(devtools)

install_github("ChaochenWang/scCDC")
```

## Usage

For detailed info on `scPNMF` method and applications, please check out the package [vignettes](https://htmlpreview.github.io/?https://github.com/ChaochenWang/scCDC/blob/main/inst/doc/scCDC.html), or with the following code in `R`: 

``` r
install_github("ChaochenWang/scCDC", build_vignettes = TRUE)
browseVignettes("scCDC")
```

## Contact

Any questions or suggestions on `scCDC` are welcomed! Please report it on [issues](https://github.com/ChaochenWang/scCDC/issues), or contact Weijian Wang (<weijian.19@intl.zju.edu.cn>), Yihui Cen (<yihui.19@intl.zju.edu.cn>) or Zezhen Lu (<zezhen.19@intl.zju.edu.cn>).

## Reference
