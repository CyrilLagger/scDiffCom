
# scDiffCom

## Overview

This package provides analysis tools to investigate changes in intercellular 
communication between two conditions of interest in scRNA-seq datasets.

It relies on a curated collection of ligand-receptor interactions (available for 
human and mouse) that have been retrieved and processed from eight public
databases.

<details>
    <summary>Click here to display the ligand-receptor sources</summary>
- [CellChat](http://www.cellchat.org/)
- [CellPhoneDB](https://www.cellphonedb.org/)
- [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/)
- [connectomeDB2020](https://github.com/forrest-lab/NATMI)
- [ICELLNET](https://github.com/soumelis-lab/ICELLNET)
- [NicheNet](https://github.com/saeyslab/nichenetr)
- [SingleCellSignalR](http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html)
- [scTensor](https://github.com/rikenbit/scTensor)
</details>&nbsp;

Using as input a [Seurat](https://satijalab.org/seurat/) object that should
contain cells annotated by cell-types and by two groups of interest
(e.g. young/old), the package infers cell-cell interactions that
potentially correspond to biological signals and significantly change between
the two conditions. For implementation details regarding the statistical
approach used in `scDiffCom` please see our (article in preparation). 


## Installation

To get the most recent version, please run:

```R
if(!require(devtools)) install.packages("devtools")
devtools::install_github("CyrilLagger/scDiffCom")
```

You can also get a stable CRAN release:

```R
#(not available yet)
```

## Typical workflow

For a usage example, please look at the vignette
available [here](https://cyrillagger.github.io/scDiffCom/)


## Reference

TBA


