
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scDiffCom

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![CRAN
status](https://www.r-pkg.org/badges/version/scDiffCom)](https://CRAN.R-project.org/package=scDiffCom)
[![Codecov test
coverage](https://codecov.io/gh/CyrilLagger/scDiffCom/branch/master/graph/badge.svg)](https://codecov.io/gh/CyrilLagger/scDiffCom?branch=master)
[![R-CMD-check](https://github.com/CyrilLagger/scDiffCom/workflows/R-CMD-check/badge.svg)](https://github.com/CyrilLagger/scDiffCom/actions)
<!-- badges: end -->

scDiffCom infers cell type to cell type communication signals from
scRNA-seq [Seurat](https://satijalab.org/seurat/) objects, and more
particularly investigates how these interactions change between two
biological conditions. The package relies on a internal collection of
ligand-receptor interactions (available for human and mouse) retrieved
from seven public and curated databases.

<details>
<summary>
Display LRI databases
</summary>

-   [CellChat](http://www.cellchat.org/)
-   [CellPhoneDB](https://www.cellphonedb.org/)
-   [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/)
-   [connectomeDB2020](https://github.com/forrest-lab/NATMI)
-   [ICELLNET](https://github.com/soumelis-lab/ICELLNET)
-   [NicheNet](https://github.com/saeyslab/nichenetr)
-   [SingleCellSignalR](http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html)

</details>

 

## Installation

``` r
# Install release version from CRAN
install.package("scDiffCom")

# Install development version from GitHub
devtools::install_github("CyrilLagger/scDiffCom")
```

## Usage

As an introduction, please look at the
[documentation](https://cyrillagger.github.io/scDiffCom/) and this
[vignette](https://cyrillagger.github.io/scDiffCom/articles/scDiffCom-vignette.html).

For a concrete and large-scale project that used scDiffCom, please look
at [scAgeCom](https://github.com/CyrilLagger/scAgeCom), our murine atlas
of age-related changes in intercellular communication.

## Citation

When using scDiffCom, please consider citing our (manuscript in
preparation).
