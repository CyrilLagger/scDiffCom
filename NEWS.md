# scDiffCom 1.1.0

* Allow to use a custom LRI database

# scDiffCom 1.0.0

01-11-2023

* Prepare for CRAN release
* Update links to the Nature Aging paper

# scDiffCom 0.2.4

11-04-2023

* Change default red-green colors in plots
* Improve README and website

# scDiffCom 0.2.3

08-05-2022

* Remove mouse and rat LRIs with "homolog_orthology_confidence" of 0, but only
if they are not already provided directly by one of the "source database"
* Remove some miss-curated LRIs

# scDiffCom 0.2.2

11-04-2022

* Compute (permutation-based) p-values and logFC for each ligand and receptor gene independently

# scDiffCom 0.2.1

20-01-2022

* Fix a bug that appeared when 'seurat_object[[]]' contained a column named 'cell_id'

# scDiffCom 0.2.0

24-11-2021

* New function `BuildShiny()` that creates a shiny app from the results
* New function `ReduceGO()` to reduce GO terms by semantic similarity

# scDiffCom 0.1.0.9000

09-11-2021

* A rat LRI database (LRI_rat) has been added (based on orthology conversion)
* `run_interaction_analysis()` accepts `rat` as species parameter
* LRI datasets have a reduced size
* Dataset `gene_ontology_level` is now available as data (instead of internal data)

# scDiffCom 0.1

04-07-2021

* All steps of `run_interaction_analysis` have been tested
* All exported functions have been tested
* LRI databases have been curated
* Complete roxygen-based documentation
* Vignette contains the main functionalities

# scDiffCom 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* Added the first proper documentation
* Added the first version of the vignette
