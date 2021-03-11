
## we expand and test some steps of the main function run_interaction_analysis ####

#the analysis can be performed in different modes according to the parameters,
#namely with or without conditions and with or without statistical test

## we create a list of different parameters corresponding to each mode ####

parameters_mode <- list(
  wrong = list(
    LRI_species = "rat",  seurat_celltype_id = c("cell_type", "cell_types"),
    seurat_condition_id = list(column_name = "age_group", cond2_name = "OLD"),
    iterations = 100.2,
    object_name = 4, seurat_assay = list("RNA"), seurat_slot = "count", log_scale = "TRUE",
    score_type = "geom", threshold_min_cells = -5, threshold_pct = 1.1,
    threshold_quantile_score = 1.5, threshold_p_value_specificity = 1.1, threshold_p_value_de = 0,
    threshold_logfc = 0, return_distributions = "FALSE", seed = 5.5, verbose = "TRUE"
  ),
  cond_stat = list(
    LRI_species = "mouse", seurat_celltype_id = "cell_type",
    seurat_condition_id = list(column_name = "age_group", cond1_name = "YOUNG", cond2_name = "OLD"),
    iterations = 10,
    object_name = "scdiffcom_cond_stat",
    seurat_assay = "RNA", seurat_slot = "data", log_scale = FALSE, score_type = "geometric_mean",
    threshold_min_cells = 5, threshold_pct = 0.1,
    threshold_quantile_score = 0.2, threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05, threshold_logfc = log(1.5),
    return_distributions = FALSE, seed = 42, verbose = FALSE
  ),
  cond_nostat = list(
    LRI_species = "mouse", seurat_celltype_id = "cell_type",
    seurat_condition_id = list(column_name = "age_group", cond1_name = "YOUNG", cond2_name = "OLD"),
    iterations = 0,
    object_name = "scdiffcom_cond_nostat",
    seurat_assay = "RNA", seurat_slot = "data", log_scale = FALSE, score_type = "geometric_mean",
    threshold_min_cells = 5, threshold_pct = 0.1,
    threshold_quantile_score = 0.2, threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05, threshold_logfc = log(1.5),
    return_distributions = FALSE, seed = 42, verbose = FALSE
  ),
  nocond_stat = list(
    LRI_species = "mouse", seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    iterations = 10,
    object_name = "scdiffcom_nocond_stat",
    seurat_assay = "RNA", seurat_slot = "data", log_scale = FALSE, score_type = "geometric_mean",
    threshold_min_cells = 5, threshold_pct = 0.1,
    threshold_quantile_score = 0.2, threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05, threshold_logfc = log(1.5),
    return_distributions = FALSE, seed = 42, verbose = FALSE
  ),
  nocond_nostat = list(
    LRI_species = "mouse", seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    iterations = 0,
    object_name = "scdiffcom_nocond_nostat",
    seurat_assay = "RNA", seurat_slot = "data", log_scale = FALSE, score_type = "geometric_mean",
    threshold_min_cells = 5, threshold_pct = 0.1,
    threshold_quantile_score = 0.2, threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05, threshold_logfc = log(1.5),
    return_distributions = FALSE, seed = 42, verbose = FALSE
  )
)

## check parameter validation #####
parameters_mode_validated <- lapply(
  parameters_mode,
  function(i) {
    validate_parameters(i, from_inputs = TRUE)$params
  }
)

parameters_mode_validated_check <- lapply(
  parameters_mode,
  function(i) {
    validate_parameters(i, from_inputs = TRUE)$check
  }
)

test_that("parameters are correclty validated", {
  lapply(parameters_mode_validated_check[2:5], expect_null)
  expect_equal(length(parameters_mode_validated$wrong), 18)
})

## small Seurat object for testing ####
seurat_test <- scDiffCom::seurat_sample_tms_liver

## check data extraction and pre-processing ####

inputs_test <- lapply(
  parameters_mode[2:5],
  function(param) {
    extract_analysis_inputs(
      seurat_object = seurat_test,
      celltype_column_id = param$seurat_celltype_id,
      sample_column_id = NULL,
      condition_column_id = param$seurat_condition_id$column_name,
      cond1_name = param$seurat_condition_id$cond1_name,
      cond2_name = param$seurat_condition_id$cond2_name,
      assay = param$seurat_assay,
      slot = param$seurat_slot,
      log_scale = param$log_scale,
      threshold_min_cells = param$threshold_min_cells,
      LRI_table = LRI_mouse$LRI_curated,
      verbose = param$verbose
    )
  }
)

test_that("data is extracted correctly in each mode", {
  expect_identical(
    inputs_test$cond_stat$data_tr["P9.MAA000907.3_11_M.1.1-1-1", "Adam15"],
    expm1(seurat_test$RNA@data["Adam15", "P9.MAA000907.3_11_M.1.1-1-1"])
  )
  expect_identical(inputs_test$cond_stat$data_tr, inputs_test$cond_nostat$data_tr)
  expect_identical(inputs_test$cond_stat$data_tr, inputs_test$nocond_stat$data_tr)
  expect_identical(inputs_test$cond_stat$data_tr, inputs_test$nocond_nostat$data_tr)
})

test_that("meta.data is extracted correctly in each mode", {
  expect_identical(
    inputs_test$cond_stat$metadata$cell_id,
    colnames(seurat_test)
  )
  expect_equivalent(
    inputs_test$cond_stat$metadata$cell_type,
    seurat_test$cell_type
  )
  expect_equivalent(
    inputs_test$cond_stat$metadata$condition,
    seurat_test$age_group
  )
  expect_identical(inputs_test$cond_stat$metadata, inputs_test$cond_nostat$metadata)
  expect_identical(inputs_test$nocond_stat$metadata, inputs_test$nocond_nostat$metadata)
})

## check the first round of the analysis ####

templates_test <- lapply(
  inputs_test,
  create_cci_template
)

cci_dt_simple_test <- lapply(
  inputs_test,
  function(inputs) {
    temp <- create_cci_template(inputs)
    run_simple_cci_analysis(
      analysis_inputs = inputs,
      cci_template = temp,
      log_scale = FALSE,
      score_type = "geometric_mean",
      threshold_min_cells = 5,
      threshold_pct = 0.1,
      compute_fast = FALSE
    )
  }
)

example_cci_test <- lapply(
  cci_dt_simple_test,
  function(i) {
    copy(i)[
      EMITTER_CELLTYPE == "hepatocyte" &
        RECEIVER_CELLTYPE == "endothelial cell of hepatic sinusoid" &
        LR_GENES == "Adam15:Itga5"]
  }
)

example_data <- list(
  EnoCond = as.vector(expm1(
    subset(
      seurat_test[c("Adam15"), ],
      subset = cell_type == "hepatocyte")$RNA@data
  )),
  RnoCond = as.vector(expm1(
    subset(
      seurat_test[c("Itga5"), ],
      subset = cell_type == "endothelial cell of hepatic sinusoid")$RNA@data
  )),
  EY = as.vector(expm1(
    subset(
      seurat_test[c("Adam15"), ],
      subset = age_group == "YOUNG" & cell_type == "hepatocyte")$RNA@data
  )),
  RY = as.vector(expm1(
    subset(
      seurat_test[c("Itga5"), ],
      subset = age_group == "YOUNG" & cell_type == "endothelial cell of hepatic sinusoid")$RNA@data
  )),
  EO = as.vector(expm1(
    subset(
      seurat_test[c("Adam15"), ],
      subset = age_group == "OLD" & cell_type == "hepatocyte")$RNA@data
  )),
  RO = as.vector(expm1(
    subset(
      seurat_test[c("Itga5"), ],
      subset = age_group == "OLD" & cell_type == "endothelial cell of hepatic sinusoid")$RNA@data
  ))
)

test_that("`simple` data.table is returned correctly", {
  expect_equivalent(mean(example_data$EnoCond > 0), example_cci_test$nocond_stat$L1_DETECTION_RATE)
  expect_equivalent(mean(example_data$RnoCond > 0), example_cci_test$nocond_stat$R1_DETECTION_RATE)
  expect_equivalent(mean(example_data$EY > 0), example_cci_test$cond_stat$L1_DETECTION_RATE_YOUNG)
  expect_equivalent(mean(example_data$RY > 0), example_cci_test$cond_stat$R1_DETECTION_RATE_YOUNG)
  expect_equivalent(mean(example_data$EO > 0), example_cci_test$cond_stat$L1_DETECTION_RATE_OLD)
  expect_equivalent(mean(example_data$RO > 0), example_cci_test$cond_stat$R1_DETECTION_RATE_OLD)
  expect_equivalent(mean(example_data$EnoCond), example_cci_test$nocond_stat$L1_EXPRESSION)
  expect_equivalent(mean(example_data$RnoCond), example_cci_test$nocond_stat$R1_EXPRESSION)
  expect_equivalent(mean(example_data$EY), example_cci_test$cond_stat$L1_EXPRESSION_YOUNG)
  expect_equivalent(mean(example_data$RY), example_cci_test$cond_stat$R1_EXPRESSION_YOUNG)
  expect_equivalent(mean(example_data$EO), example_cci_test$cond_stat$L1_EXPRESSION_OLD)
  expect_equivalent(mean(example_data$RO), example_cci_test$cond_stat$R1_EXPRESSION_OLD)
  expect_equivalent(sqrt(mean(example_data$EnoCond)*mean(example_data$RnoCond)), example_cci_test$nocond_stat$CCI_SCORE)
  expect_equivalent(sqrt(mean(example_data$EnoCond)*mean(example_data$RnoCond)), example_cci_test$nocond_nostat$CCI_SCORE)
  expect_equivalent(sqrt(mean(example_data$EY)*mean(example_data$RY)), example_cci_test$cond_stat$CCI_SCORE_YOUNG)
  expect_equivalent(sqrt(mean(example_data$EO)*mean(example_data$RO)), example_cci_test$cond_stat$CCI_SCORE_OLD)
  expect_identical(cci_dt_simple_test$cond_stat, cci_dt_simple_test$cond_nostat)
  expect_identical(cci_dt_simple_test$nocond_stat, cci_dt_simple_test$nocond_nostat)
})

## check overall permutation analysis ####


cci_permutation_test <- lapply(
  c(cond = 1, nocond = 3),
  function(i) {
    set.seed(42)
    run_stat_analysis(
      analysis_inputs = inputs_test[[i]],
      cci_dt_simple = cci_dt_simple_test[[i]],
      iterations = 10,
      return_distributions = FALSE,
      score_type = "geometric_mean",
      verbose = FALSE
    )
  }
)

cci_permutation_test_with_distr <- lapply(
  c(cond = 1, nocond = 3),
  function(i) {
    set.seed(42)
    run_stat_analysis(
      analysis_inputs = inputs_test[[i]],
      cci_dt_simple = cci_dt_simple_test[[i]],
      iterations = 10,
      return_distributions = TRUE,
      score_type = "geometric_mean",
      verbose = FALSE
    )
  }
)

test_that("data.table returned by permutation test is correctly formatted", {
  expect_identical(cci_dt_simple_test$cond_stat$LR_GENES, cci_permutation_test$cond$cci_raw$LR_GENES)
  expect_identical(cci_dt_simple_test$cond_stat[, 4:38], cci_permutation_test$cond$cci_raw[, 4:38])
  expect_identical(cci_dt_simple_test$nocond_stat[, 4:22], cci_permutation_test$nocond$cci_raw[, 4:22])
  expect_identical(cci_dt_simple_test$cond_stat$LR_GENES, cci_permutation_test_with_distr$cond$cci_raw$LR_GENES)
  expect_identical(cci_dt_simple_test$cond_stat[, 4:38], cci_permutation_test_with_distr$cond$cci_raw[, 4:38])
  expect_identical(cci_dt_simple_test$nocond_stat[, 4:22], cci_permutation_test_with_distr$nocond$cci_raw[, 4:22])
})

test_that("permutation results are identical (using same seed) with or without returning the disributions", {
  expect_identical(cci_permutation_test_with_distr$cond$cci_raw$P_VALUE_DE,
                   cci_permutation_test$cond$cci_raw$P_VALUE_DE )
  expect_identical(cci_permutation_test_with_distr$cond$cci_raw$P_VALUE_OLD,
                   cci_permutation_test$cond$cci_raw$P_VALUE_OLD)
  expect_identical(cci_permutation_test_with_distr$cond$cci_raw$P_VALUE_YOUNG,
                   cci_permutation_test$cond$cci_raw$P_VALUE_YOUNG)
})


## benchmarking of permutation test and internal functions ####
#library(profvis)
#library(microbenchmark)
#toDO

## check object is returned correctly ####

scdiffcom_objects <- lapply(
  parameters_mode_validated[2:5],
  function(params) {
    run_internal_raw_analysis(
      seurat_object = seurat_test,
      LRI_table = LRI_mouse$LRI_curated,
      params = params
    )
  }
)

test_that("scDiffCom raw object is returned properly", {
  lapply(scdiffcom_objects,
         function(object) {expect_s4_class(object, "scDiffCom")}
  )
})

## Check filtering without ORA ####

scdiffcom_objects <- lapply(
  scdiffcom_objects,
  function(object) {
    run_filtering_and_ora(
      object  = object,
      new_threshold_quantile_score = NULL,
      new_threshold_p_value_specificity = NULL,
      new_threshold_p_value_de = NULL,
      new_threshold_logfc = NULL,
      skip_ora = TRUE,
      verbose = FALSE,
      class_signature = "scDiffCom"
    )
  }
)

#test_that("todo", {
#
#})


## check ORA ####

scdiffcom_objects <- lapply(
  scdiffcom_objects,
  function(object) {
    run_ora(
      object  = object,
      categories = c("ER_CELLTYPES", "LR_GENES", "GO_TERMS", "KEGG_PWS"),
      overwrite = TRUE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      verbose = TRUE,
      class_signature = "scDiffCom",
      global = FALSE
    )
  }
)

# ORA on LR_GENES

contingency_table_test_LR_up <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES == "Apoe:Ldlr" & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES == "Apoe:Ldlr" & REGULATION != "UP", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES != "Apoe:Ldlr" & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES != "Apoe:Ldlr" & REGULATION != "UP", .N]
  ),
  2,2
)

contingency_table_test_LR_down <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES == "Csf2:Itgb1" & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES == "Csf2:Itgb1" & REGULATION != "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES != "Csf2:Itgb1" & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES != "Csf2:Itgb1" & REGULATION != "DOWN", .N]
  ),
  2,2
)

contingency_table_test_LR_flat <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES == "Adam10:Axl" & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES == "Adam10:Axl" & REGULATION != "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES != "Adam10:Axl" & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES != "Adam10:Axl" & REGULATION != "FLAT", .N]
  ),
  2,2
)

fisher_LR_up <- fisher.test(contingency_table_test_LR_up)
fisher_LR_down <- fisher.test(contingency_table_test_LR_down)
fisher_LR_flat <- fisher.test(contingency_table_test_LR_flat)


test_that("fisher test is done correctly on LRIs", {
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$LR_GENES[VALUE == "Apoe:Ldlr"]$OR_UP,
    fisher_LR_up$estimate
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$LR_GENES[VALUE == "Apoe:Ldlr"]$P_VALUE_UP,
    fisher_LR_up$p.value
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$LR_GENES[VALUE == "Csf2:Itgb1"]$OR_DOWN,
    fisher_LR_down$estimate
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$LR_GENES[VALUE == "Csf2:Itgb1"]$P_VALUE_DOWN,
    fisher_LR_down$p.value
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$LR_GENES[VALUE == "Adam10:Axl"]$OR_FLAT,
    fisher_LR_flat$estimate
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$LR_GENES[VALUE == "Adam10:Axl"]$P_VALUE_FLAT,
    fisher_LR_flat$p.value
  )
})

#ORA on GO terms

genes_GO_UP <- LRI_mouse$LRI_curated_GO[GO_ID == "GO:0002376"]$LR_GENES
contingency_table_test_GO_up <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES %in% genes_GO_UP & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES %in% genes_GO_UP & REGULATION != "UP", .N],
    scdiffcom_objects$cond_stat@cci_detected[!(LR_GENES %in% genes_GO_UP) & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_detected[!(LR_GENES %in% genes_GO_UP) & REGULATION != "UP", .N]
  ),
  2,2
)
fisher_GO_up <- fisher.test(contingency_table_test_GO_up)

genes_GO_DOWN <- LRI_mouse$LRI_curated_GO[GO_ID == "GO:0031625"]$LR_GENES
contingency_table_test_GO_down <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES %in% genes_GO_DOWN & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES %in% genes_GO_DOWN & REGULATION != "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_detected[!(LR_GENES %in% genes_GO_DOWN) & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_detected[!(LR_GENES %in% genes_GO_DOWN) & REGULATION != "DOWN", .N]
  ),
  2,2
)
fisher_GO_down <- fisher.test(contingency_table_test_GO_down)

genes_GO_FLAT <- LRI_mouse$LRI_curated_GO[GO_ID == "GO:0022406"]$LR_GENES
contingency_table_test_GO_flat <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES %in% genes_GO_FLAT & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_detected[LR_GENES %in% genes_GO_FLAT & REGULATION != "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_detected[!(LR_GENES %in% genes_GO_FLAT) & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_detected[!(LR_GENES %in% genes_GO_FLAT) & REGULATION != "FLAT", .N]
  ),
  2,2
)
fisher_GO_flat <- fisher.test(contingency_table_test_GO_flat)
fisher_GO_flat_cor <- fisher.test(contingency_table_test_GO_flat+1)

test_that("fisher test is done correctly on GO terms", {
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$GO_TERMS[VALUE_BIS == "GO:0002376"]$OR_UP,
    fisher_GO_up$estimate
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$GO_TERMS[VALUE_BIS == "GO:0002376"]$P_VALUE_UP,
    fisher_GO_up$p.value
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$GO_TERMS[VALUE_BIS == "GO:0031625"]$OR_DOWN,
    fisher_GO_down$estimate
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$GO_TERMS[VALUE_BIS == "GO:0031625"]$P_VALUE_DOWN,
    fisher_GO_down$p.value
  )
  expect_equal(
    scdiffcom_objects$cond_stat@ora_default$GO_TERMS[VALUE_BIS == "GO:0022406"]$OR_FLAT[[1]],
    fisher_GO_flat_cor$estimate[[1]],
    tolerance = 0.01
  )
  expect_equivalent(
    scdiffcom_objects$cond_stat@ora_default$GO_TERMS[VALUE_BIS == "GO:0022406"]$P_VALUE_FLAT,
    fisher_GO_flat$p.value
  )
})

## Check BuildNetwork step by step ####

## Check construct_graph ####

# ig_cond_stat <- construct_graph(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   ora_default_ER = scdiffcom_objects$cond_stat@ora_default$ER_CELLTYPES
# )
# plot(ig_cond_stat)
# igraph::V(ig_cond_stat)
# igraph::E(ig_cond_stat)
# igraph::edge.attributes(ig_cond_stat)
# igraph::vertex.attributes(ig_cond_stat)
#
# ig_nocond_stat <- construct_graph(
#   cci_detected = scdiffcom_objects$nocond_stat@cci_detected,
#   conds = NULL,
#   ora_default_ER = NULL
# )
# plot(ig_nocond_stat)
# igraph::V(ig_nocond_stat)
# igraph::E(ig_nocond_stat)
# igraph::edge.attributes(ig_nocond_stat)
# igraph::vertex.attributes(ig_nocond_stat)
#
# ## Check setup_graph #####
#
# setup_test <- setup_graph(
#   G = ig_cond_stat,
#   conds = c("YOUNG", "OLD"),
#   config = setup_graph_config() ,
#   disperse = TRUE
# )
#
#
# ## Check build_igraph ####
#
# igb_cond_stat <- build_igraph(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   ora_default_ER = scdiffcom_objects$cond_stat@ora_default$ER_CELLTYPES,
#   network_layout_type = "bipartite"
# )
# plot(igb_cond_stat)
#
# ## Check network_skeleton ####
#
# interactive_from_igraph <- function(
#   cci_detected,
#   ora_default_LR,
#   object_name,
#   G,
#   network_representation_type,
#   network_layout_type
# )
#
# test_nodes <- setDT(igraph::as_data_frame(igb_cond_stat, what = "vertices"))
# test_edges <- setDT(igraph::as_data_frame(igb_cond_stat, what = "edges"))
# test_edges <- test_edges[ORA_TYPE != "NONE"]
#
# test <- build_network_skeleton(
#   test_nodes,
#   test_edges,
#   "ORA",
#   "bipartite",
#   setup_graph_config(),
#   "test"
# )
# test
#
#
# types_of_network <- c(
#     "condition1_network",
#     "condition2_network",
#     "difference_network",
#     "up_regulated_network",
#     "down_regulated_network",
#     "ORA_network"
# )
#
# layouts_of_network = c(
#   "conventional",
#   "bipartite"
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "ORA_network",
#   layout_type = "bipartite",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "ORA_network",
#   layout_type = "conventional",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "difference_network",
#   layout_type = "bipartite",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "difference_network",
#   layout_type = "conventional",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "down_regulated_network",
#   layout_type = "bipartite",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "condition2_network",
#   layout_type = "bipartite",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "condition2_network",
#   layout_type = "conventional",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$nocond_stat,
#   network_type = "condition1_network",
#   layout_type = "bipartite",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$nocond_stat,
#   network_type = "condition1_network",
#   layout_type = "conventional",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# build_interactive_network(
#   object = scdiffcom_objects$cond_stat,
#   network_type = "up_regulated_network",
#   layout_type = "conventional",
#   class_signature = "scdiffCom",
#   subobject_name = NULL
# )
#
# edge_test <- build_edge_table(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   ora_default_ER = scdiffcom_objects$cond_stat@ora_default$ER_CELLTYPES,
#   ora_default_LR = scdiffcom_objects$cond_stat@ora_default$LR_GENES,
#   network_representation_type = "ORA",
#   network_layout_type = "bipartite",
#   config = setup_graph_config()
# )
#
# vertex_test <- build_vertex_table(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   layout_type = "bipartite",
#   config = setup_graph_config()
# )
# vertex_test
#
#
#
# test <- extract_vertex_metadata(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   layout_type = "bipartite"
# )
#
# test
#
# #
# test <- build_igraph(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   ora_default_ER = scdiffcom_objects$cond_stat@ora_default$ER_CELLTYPES,
#   ora_default_LR = scdiffcom_objects$cond_stat@ora_default$LR_GENES,
#   network_type = "ORA_network",
#   layout_type = "bipartite",
#   config = setup_graph_config()
# )
# plot(test)
#
# n_emit <- igraph::V(test)$name[igraph::V(test)$vertex_types]
# edg <- setDT(igraph::edge.attributes(test))
# edg

# test
# plot(test, edge.loop.angle = igraph::E(test)$edge.loop.angle )
#
#
# test <- interactive_from_igraph(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   ora_default_ER = scdiffcom_objects$cond_stat@ora_default$ER_CELLTYPES,
#   ora_default_LR = scdiffcom_objects$cond_stat@ora_default$LR_GENES,
#   network_representation_type = "ORA",
#   network_layout_type = "bipartite",
#   object_name = "test"
# )
#
# edge_test <- build_edge_table(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD"),
#   ora_default_ER = scdiffcom_objects$cond_stat@ora_default$ER_CELLTYPES,
#   ora_default_LR = scdiffcom_objects$cond_stat@ora_default$LR_GENES,
#   network_type = "ORA_network",
#   layout_type = "bipartite",
#   config = setup_graph_config()
# )
#
# sort_bipartite_vertices(
#   test,
#   edge_test,
#   "ORA_network"
# )

#
# test <- extract_edge_metadata(
#   cci_detected = scdiffcom_objects$cond_stat@cci_detected,
#   conds = c("YOUNG", "OLD")
# )
# test <- process_celltype_pairs_enrichment(
#   ora_ER_cells =
# )

# scdiffcom_objects$cond_stat@parameters$object_name <- "abc"
#
# scdiffcom_test_comb <- Combine_scDiffCom(
#   l = list(scdiffcom_objects$cond_stat, scdiffcom_objects2$cond_stat),
#   object_name = "test",
#   verbose = TRUE
# )
#
# BuildNetwork(
#   scdiffcom_test_comb,
#   ID = "abc",
#   abbreviation_table = NULL
# )
