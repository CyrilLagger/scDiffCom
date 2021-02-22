
## we expand and test some steps of the main function run_interaction_analysis ####

#the analysis can be performed in different modes according to the parameters,
#namely with or without conditions and with or without statistical test

## we create a list of different parameters corresponding to each mode ####

parameters_mode <- list(
  wrong = list(
    LRdb_species = "rat",  seurat_celltype_id = c("cell_type", "cell_types"),
    seurat_condition_id = list(column_name = "age_group", cond2_name = "OLD"),
    iterations = 100.2,
    object_name = 4, seurat_assay = list("RNA"), seurat_slot = "count", log_scale = "TRUE",
    score_type = "geom", threshold_min_cells = -5, threshold_pct = 1.1,
    threshold_quantile_score = 1.5, threshold_p_value_specificity = 1.1, threshold_p_value_de = 0,
    threshold_logfc = 0, return_distributions = "FALSE", seed = 5.5, verbose = "TRUE"
  ),
  cond_stat = list(
    LRdb_species = "mouse", seurat_celltype_id = "cell_type",
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
    LRdb_species = "mouse", seurat_celltype_id = "cell_type",
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
    LRdb_species = "mouse", seurat_celltype_id = "cell_type",
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
    LRdb_species = "mouse", seurat_celltype_id = "cell_type",
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

## check data extraction and preprocessing ####

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
      LRdb_table = LRdb_mouse$LRdb_curated,
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

example_cci_test <- copy(cci_dt_simple_test$cond_stat)[
  EMITTER_CELLTYPE == "hepatocyte" &
  RECEIVER_CELLTYPE == "endothelial cell of hepatic sinusoid" &
    LR_GENES == "Adam15:Itga5"]

example_data_EY <- as.vector(expm1(
  subset(
    seurat_test[c("Adam15"), ],
    subset = age_group == "YOUNG" & cell_type == "hepatocyte")$RNA@data
))
example_data_RY <- as.vector(expm1(
  subset(
    seurat_test[c("Itga5"), ],
    subset = age_group == "YOUNG" & cell_type == "endothelial cell of hepatic sinusoid")$RNA@data
))
example_data_EO <- as.vector(expm1(
  subset(
    seurat_test[c("Adam15"), ],
    subset = age_group == "OLD" & cell_type == "hepatocyte")$RNA@data
))
example_data_RO <- as.vector(expm1(
  subset(
    seurat_test[c("Itga5"), ],
    subset = age_group == "OLD" & cell_type == "endothelial cell of hepatic sinusoid")$RNA@data
))

test_that("`simple` data.table is returned correctly", {
  expect_equivalent(mean(example_data_EY > 0), example_cci_test$L1_DETECTION_RATE_YOUNG)
  expect_equivalent(mean(example_data_RY > 0), example_cci_test$R1_DETECTION_RATE_YOUNG)
  expect_equivalent(mean(example_data_EO > 0), example_cci_test$L1_DETECTION_RATE_OLD)
  expect_equivalent(mean(example_data_RO > 0), example_cci_test$R1_DETECTION_RATE_OLD)
  expect_equivalent(mean(example_data_EY), example_cci_test$L1_EXPRESSION_YOUNG)
  expect_equivalent(mean(example_data_RY), example_cci_test$R1_EXPRESSION_YOUNG)
  expect_equivalent(mean(example_data_EO), example_cci_test$L1_EXPRESSION_OLD)
  expect_equivalent(mean(example_data_RO), example_cci_test$R1_EXPRESSION_OLD)
  expect_equivalent(sqrt(mean(example_data_EY)*mean(example_data_RY)), example_cci_test$CCI_SCORE_YOUNG)
  expect_equivalent(sqrt(mean(example_data_EO)*mean(example_data_RO)), example_cci_test$CCI_SCORE_OLD)
  expect_identical(cci_dt_simple_test$cond_stat, cci_dt_simple_test$cond_nostat)
  expect_identical(cci_dt_simple_test$nocond_stat, cci_dt_simple_test$nocond_nostat)
})

## check the permutation analysis ####
cci_permutation_test <- lapply(
  c(cond = 1, nocond = 3),
  function(i) {
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

test_that("permutation test is performed correclty", {
  expect_identical(cci_dt_simple_test$cond_stat$LR_GENES, cci_permutation_test$cond$cci_raw$LR_GENES)
  expect_identical(cci_dt_simple_test$cond_stat[, 4:38], cci_permutation_test$cond$cci_raw[, 4:38])
  expect_identical(cci_dt_simple_test$nocond_stat[, 4:22], cci_permutation_test$nocond$cci_raw[, 4:22])
})

## check object is returned correclty ####

scdiffcom_objects <- lapply(
  parameters_mode_validated[2:5],
  function(params) {
    run_internal_raw_analysis(
      seurat_object = seurat_test,
      LRdb_table = LRdb_mouse$LRdb_curated,
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
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      verbose = TRUE,
      class_signature = "scDiffCom",
      global = FALSE
    )
  }
)

scdiffcom_objects$cond_stat@ora_default$LR_GENES





