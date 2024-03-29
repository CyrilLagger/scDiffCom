
data.table::setDTthreads(2)

## We test key functions of the package ####
# roughly in the order how they are called in a standard analysis after
# calling run_interaction_analysis

## Create a list of different parameters corresponding to each mode ####

parameters_mode <- list(
  wrong = list(
    LRI_species = "cat",
    seurat_celltype_id = c("cell_type", "cell_types"),
    seurat_condition_id = list(
      column_name = "age_group",
      cond2_name = "OLD"
    ),
    iterations = 100.2,
    object_name = 4,
    seurat_assay = list("RNA"),
    seurat_slot = "count",
    log_scale = "TRUE",
    score_type = "geom",
    threshold_min_cells = -5,
    threshold_pct = 1.1,
    threshold_quantile_score = 1.5,
    threshold_p_value_specificity = 1.1,
    threshold_p_value_de = 0,
    threshold_logfc = 0,
    return_distributions = "FALSE",
    seed = 5.5,
    verbose = "TRUE",
    custom_LRI_tables = "something_wrong"
  ),
  cond_stat = list(
    LRI_species = "mouse",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "age_group",
      cond1_name = "YOUNG",
      cond2_name = "OLD"
    ),
    iterations = 10,
    object_name = "scdiffcom_cond_stat",
    seurat_assay = "RNA",
    seurat_slot = "data",
    log_scale = FALSE,
    score_type = "geometric_mean",
    threshold_min_cells = 5,
    threshold_pct = 0.1,
    threshold_quantile_score = 0.2,
    threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05,
    threshold_logfc = log(1.5),
    return_distributions = FALSE,
    seed = 42,
    verbose = TRUE,
    custom_LRI_tables = NULL
  ),
  cond_nostat = list(
    LRI_species = "mouse",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = list(
      column_name = "age_group",
      cond1_name = "YOUNG",
      cond2_name = "OLD"
    ),
    iterations = 0,
    object_name = "scdiffcom_cond_nostat",
    seurat_assay = "RNA",
    seurat_slot = "data",
    log_scale = FALSE,
    score_type = "geometric_mean",
    threshold_min_cells = 5,
    threshold_pct = 0.1,
    threshold_quantile_score = 0.2,
    threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05,
    threshold_logfc = log(1.5),
    return_distributions = FALSE,
    seed = 42,
    verbose = FALSE,
    custom_LRI_tables = NULL
  ),
  nocond_stat = list(
    LRI_species = "mouse",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    iterations = 10,
    object_name = "scdiffcom_nocond_stat",
    seurat_assay = "RNA",
    seurat_slot = "data",
    log_scale = FALSE,
    score_type = "geometric_mean",
    threshold_min_cells = 5,
    threshold_pct = 0.1,
    threshold_quantile_score = 0.2,
    threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05,
    threshold_logfc = log(1.5),
    return_distributions = FALSE,
    seed = 42,
    verbose = FALSE,
    custom_LRI_tables = NULL
  ),
  nocond_nostat = list(
    LRI_species = "mouse",
    seurat_celltype_id = "cell_type",
    seurat_condition_id = NULL,
    iterations = 0,
    object_name = "scdiffcom_nocond_nostat",
    seurat_assay = "RNA",
    seurat_slot = "data",
    log_scale = FALSE,
    score_type = "geometric_mean",
    threshold_min_cells = 5,
    threshold_pct = 0.1,
    threshold_quantile_score = 0.2,
    threshold_p_value_specificity = 0.05,
    threshold_p_value_de = 0.05,
    threshold_logfc = log(1.5),
    return_distributions = FALSE,
    seed = 42,
    verbose = FALSE,
    custom_LRI_tables = NULL
  )
)

## Check input parameter validation #####

parameters_mode_validated <- lapply(
  parameters_mode,
  function(i) {
    validate_parameters(
      i,
      from_inputs = TRUE
    )$params
  }
)

parameters_mode_validated_check <- lapply(
  parameters_mode,
  function(i) {
    validate_parameters(
      i,
      from_inputs = TRUE
    )$check
  }
)

test_that(
  "parameters are correclty validated",
  {
    lapply(
      parameters_mode_validated_check[2:5],
      expect_null
    )
    expect_length(
      parameters_mode_validated$wrong,
      19
    )
  }
)

## Load a light Seurat object for testing ####

seurat_test <- scDiffCom::seurat_sample_tms_liver

test_that(
  "sample object is of class Seurat",
  {
    expect_s4_class(seurat_test, "Seurat")
  }
)

## Check data extraction and pre-processing ####

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
      LRI_species = "mouse",
      verbose = param$verbose
    )
  }
)

test_that(
  "data is extracted correctly in each mode",
  {
    expect_identical(
      inputs_test$cond_stat$data_tr["P9.MAA000907.3_11_M.1.1-1-1", "Adam15"],
      expm1(seurat_test$RNA@data["Adam15", "P9.MAA000907.3_11_M.1.1-1-1"])
    )
    expect_identical(
      inputs_test$cond_stat$data_tr,
      inputs_test$cond_nostat$data_tr
    )
    expect_identical(
      inputs_test$cond_stat$data_tr,
      inputs_test$nocond_stat$data_tr
    )
    expect_identical(
      inputs_test$cond_stat$data_tr,
      inputs_test$nocond_nostat$data_tr
    )
  }
)

test_that(
  "meta.data is extracted correctly in each mode",
  {
    expect_identical(
      inputs_test$cond_stat$metadata$cell_id,
      colnames(seurat_test)
    )
    expect_equal(
      inputs_test$cond_stat$metadata$cell_type,
      seurat_test$cell_type,
      ignore_attr = TRUE
    )
    expect_equal(
      inputs_test$cond_stat$metadata$condition,
      seurat_test$age_group,
      ignore_attr = TRUE
    )
    expect_identical(
      inputs_test$cond_stat$metadata,
      inputs_test$cond_nostat$metadata
    )
    expect_identical(
      inputs_test$nocond_stat$metadata,
      inputs_test$nocond_nostat$metadata
    )
  }
)

test_that(
  "cell types are extracted correctly in each mode",
  {
    lapply(
      inputs_test,
      function(i) {
        expect_identical(
          sort(i$cell_types),
          sort(unique(seurat_test$cell_type))
        )
      }
    )
  }
)

test_that(
  "condition is extracted correctly in the relevant mode",
  {
    expect_identical(
      inputs_test$cond_stat$condition,
      list(is_cond = TRUE, is_samp = FALSE, cond1 = "YOUNG", cond2 = "OLD")
    )
  }
)

## Check the first round of the analysis ####

templates_test <- lapply(
  inputs_test,
  create_cci_template
)

test_that(
  "cci_template returns the correct number of CCIs",
  {
    lapply(
      templates_test,
      function(i) {
        expect_equal(
          nrow(i),
          29325
        )
      }
    )
  }
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
        LRI == "Adam15:Itga5"
    ]
  }
)

example_data <- list(
  EnoCond = as.vector(
    expm1(
      subset(
        seurat_test[c("Adam15"), ],
        subset = cell_type == "hepatocyte"
      )$RNA@data
    )
  ),
  RnoCond = as.vector(
    expm1(
      subset(
        seurat_test[c("Itga5"), ],
        subset = cell_type == "endothelial cell of hepatic sinusoid"
      )$RNA@data
    )
  ),
  EY = as.vector(
    expm1(
      subset(
        seurat_test[c("Adam15"), ],
        subset = age_group == "YOUNG" & cell_type == "hepatocyte"
      )$RNA@data
    )
  ),
  RY = as.vector(
    expm1(
      subset(
        seurat_test[c("Itga5"), ],
        subset = age_group == "YOUNG" &
          cell_type == "endothelial cell of hepatic sinusoid"
      )$RNA@data
    )
  ),
  EO = as.vector(
    expm1(
      subset(
        seurat_test[c("Adam15"), ],
        subset = age_group == "OLD" &
          cell_type == "hepatocyte"
      )$RNA@data
    )
  ),
  RO = as.vector(
    expm1(
      subset(
        seurat_test[c("Itga5"), ],
        subset = age_group == "OLD" &
          cell_type == "endothelial cell of hepatic sinusoid"
      )$RNA@data
    )
  )
)

test_that(
  "`simple` data.table is returned correctly",
  {
    expect_equal(
      mean(
        example_data$EnoCond > 0
      ),
      example_cci_test$nocond_stat$L1_DETECTION_RATE,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$RnoCond > 0
      ),
      example_cci_test$nocond_stat$R1_DETECTION_RATE,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$EY > 0
      ),
      example_cci_test$cond_stat$L1_DETECTION_RATE_YOUNG,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$RY > 0
      ),
      example_cci_test$cond_stat$R1_DETECTION_RATE_YOUNG,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$EO > 0
      ),
      example_cci_test$cond_stat$L1_DETECTION_RATE_OLD,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$RO > 0
      ),
      example_cci_test$cond_stat$R1_DETECTION_RATE_OLD,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$EnoCond
      ),
      example_cci_test$nocond_stat$L1_EXPRESSION,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$RnoCond
        ),
      example_cci_test$nocond_stat$R1_EXPRESSION,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$EY
      ),
      example_cci_test$cond_stat$L1_EXPRESSION_YOUNG,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$RY
      ),
      example_cci_test$cond_stat$R1_EXPRESSION_YOUNG,
      ignore_attr = TRUE
    )
    expect_equal(
      mean(
        example_data$EO
        ), example_cci_test$cond_stat$L1_EXPRESSION_OLD,
      ignore_attr = TRUE
      )
    expect_equal(
      mean(
        example_data$RO
      ),
      example_cci_test$cond_stat$R1_EXPRESSION_OLD,
      ignore_attr = TRUE
    )
    expect_equal(
      sqrt(
        mean(example_data$EnoCond)*mean(example_data$RnoCond)
      ),
      example_cci_test$nocond_stat$CCI_SCORE,
      ignore_attr = TRUE
    )
    expect_equal(
      sqrt(
        mean(example_data$EnoCond)*mean(example_data$RnoCond)
      ),
      example_cci_test$nocond_nostat$CCI_SCORE,
      ignore_attr = TRUE
    )
    expect_equal(
      sqrt(
        mean(example_data$EY)*mean(example_data$RY)
      ),
      example_cci_test$cond_stat$CCI_SCORE_YOUNG,
      ignore_attr = TRUE
    )
    expect_equal(
      sqrt(
        mean(example_data$EO)*mean(example_data$RO)
      ),
      example_cci_test$cond_stat$CCI_SCORE_OLD,
      ignore_attr = TRUE
    )
    expect_identical(
      cci_dt_simple_test$cond_stat,
      cci_dt_simple_test$cond_nostat
    )
    expect_identical(
      cci_dt_simple_test$nocond_stat,
      cci_dt_simple_test$nocond_nostat
    )
  }
)

## Check permutation analysis ####

cci_permutation_test <- lapply(
  c(cond = 1, nocond = 3),
  function(i) {
    set.seed(123)
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
    set.seed(123)
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

test_that(
  "data.table returned by permutation test is correctly formatted", {
    expect_identical(
      cci_dt_simple_test$cond_stat$LRI,
      cci_permutation_test$cond$cci_raw$LRI
    )
    expect_identical(
      cci_dt_simple_test$cond_stat[, 4:38],
      cci_permutation_test$cond$cci_raw[, 4:38]
    )
    expect_identical(
      cci_dt_simple_test$nocond_stat[, 4:22],
      cci_permutation_test$nocond$cci_raw[, 4:22]
    )
    expect_identical(
      cci_dt_simple_test$cond_stat$LRI,
      cci_permutation_test_with_distr$cond$cci_raw$LRI
    )
    expect_identical(
      cci_dt_simple_test$cond_stat[, 4:38],
      cci_permutation_test_with_distr$cond$cci_raw[, 4:38]
    )
    expect_identical(
      cci_dt_simple_test$nocond_stat[, 4:22],
      cci_permutation_test_with_distr$nocond$cci_raw[, 4:22]
    )
  }
)

test_that(
  "permutation results are identical (using same seed) with or without returning the distributions", {
    expect_identical(
      cci_permutation_test_with_distr$cond$cci_raw$P_VALUE_DE,
      cci_permutation_test$cond$cci_raw$P_VALUE_DE
    )
    expect_identical(
      cci_permutation_test_with_distr$cond$cci_raw$P_VALUE_OLD,
      cci_permutation_test$cond$cci_raw$P_VALUE_OLD
    )
    expect_identical(
      cci_permutation_test_with_distr$cond$cci_raw$P_VALUE_YOUNG,
      cci_permutation_test$cond$cci_raw$P_VALUE_YOUNG
    )
  }
)

## Check object is returned correctly ####

scdiffcom_objects <- lapply(
  parameters_mode_validated[2:5],
  function(params) {
    run_internal_raw_analysis(
      seurat_object = seurat_test,
      LRI_table = LRI_mouse$LRI_curated,
      LRI_species = "mouse",
      params = params
    )
  }
)

test_that(
  "scDiffCom raw object is returned properly",
  {
    lapply(
      scdiffcom_objects,
      function(object) {
        expect_s4_class(object, "scDiffCom")
      }
    )
  }
)

## Check main function is not running if not working on Seurat object ####

test_that(
  "run_interaction_analysis returns error if not passed a Seurat object", {
    expect_error(
      run_interaction_analysis(
        seurat_object = "not a Seurat object",
        LRI_species = "mouse",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        )
      )
    )
  }
)

## Check main function is running ####

test_that(
  "run_interaction_analysis is running correctly", {
    expect_s4_class(
      run_interaction_analysis(
        seurat_object = seurat_test,
        LRI_species = "mouse",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        )
      ),
     "scDiffCom"
    )
  }
)

## Check the use of custom_LRI_table ####

test_that(
  "run_interaction_analysis returns error if custom_LRI_tables$LRI is not a data.table", {
    expect_error(
      run_interaction_analysis(
        seurat_object = seurat_test,
        LRI_species = "mouse",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        ),
        custom_LRI_tables = "something_wrong"
      )
    )
  }
)

test_that(
  "run_interaction_analysis returns error if custom_LRI_tables is not a data.table", {
    expect_error(
      run_interaction_analysis(
        seurat_object = seurat_test,
        LRI_species = "mouse",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        ),
        custom_LRI_tables = list(LRI = data.frame(x = 1:4))
      )
    )
  }
)

test_that(
  "run_interaction_analysis returns error if custom_LRI_table is not formatted correctly", {
    expect_error(
      run_interaction_analysis(
        seurat_object = seurat_test,
        LRI_species = "custom",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        ),
        custom_LRI_tables = list(LRI = data.table(x = 1:4))
      )
    )
  }
)

test_that(
  "error if LRI_species and custom_LRI_table are inconsistent", {
    expect_error(
      run_interaction_analysis(
        seurat_object = seurat_test,
        LRI_species = "mouse",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        ),
        custom_LRI_table = list(LRI = LRI_mouse$LRI_curated[1:500, ])
      )
    )
  }
)

test_that(
  "run_interaction_analysis is running correctly with custom LRIs and no annotations", {
    expect_s4_class(
      run_interaction_analysis(
        seurat_object = seurat_test,
        LRI_species = "custom",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        ),
        custom_LRI_table = list(LRI = LRI_mouse$LRI_curated[1:500, ])
      ),
     "scDiffCom"
    )
  }
)

custom_LRI <- LRI_mouse$LRI_curated[1:500, ]
custom_LRI_GO <- LRI_mouse$LRI_curated_GO[LRI %in% custom_LRI$LRI]

test_that(
  "run_interaction_analysis is running correctly with custom LRIs and annotations", {
    expect_s4_class(
      run_interaction_analysis(
        seurat_object = seurat_test,
        LRI_species = "custom",
        seurat_celltype_id = "cell_type",
        seurat_condition_id = list(
          column_name = "age_group",
          cond1_name = "YOUNG",
          cond2_name = "OLD"
        ),
        custom_LRI_table = list(LRI = custom_LRI, my_GO = custom_LRI_GO)
      ),
     "scDiffCom"
    )
  }
)

## Check internal filtering without ORA ####

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

# nothing particular to test as long as the previous function runs

## Check internal ORA ####

scdiffcom_objects <- lapply(
  scdiffcom_objects,
  function(object) {
    run_ora(
      object  = object,
      categories = c(
        "LRI",
        "LIGAND_COMPLEX",
        "RECEPTOR_COMPLEX",
        "ER_CELLTYPES",
        "EMITTER_CELLTYPE",
        "RECEIVER_CELLTYPE",
        "GO_TERMS",
        "KEGG_PWS"
      ),
      extra_annotations = NULL,
      overwrite = TRUE,
      stringent_or_default = "default",
      stringent_logfc_threshold = NULL,
      verbose = TRUE,
      class_signature = "scDiffCom",
      global = FALSE
    )
  }
)

# ORA double-check on LRIs

contingency_table_test_LR_up <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_table_detected[LRI == "Apoe:Ldlr" & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI == "Apoe:Ldlr" & REGULATION != "UP", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI != "Apoe:Ldlr" & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI != "Apoe:Ldlr" & REGULATION != "UP", .N]
  ),
  2,
  2
)

contingency_table_test_LR_down <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_table_detected[LRI == "Csf2:Itgb1" & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI == "Csf2:Itgb1" & REGULATION != "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI != "Csf2:Itgb1" & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI != "Csf2:Itgb1" & REGULATION != "DOWN", .N]
  ),
  2,
  2
)

contingency_table_test_LR_flat <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_table_detected[LRI == "Adam10:Axl" & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI == "Adam10:Axl" & REGULATION != "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI != "Adam10:Axl" & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI != "Adam10:Axl" & REGULATION != "FLAT", .N]
  ),
  2,
  2
)

fisher_LR_up <- fisher.test(contingency_table_test_LR_up)
fisher_LR_down <- fisher.test(contingency_table_test_LR_down)
fisher_LR_flat <- fisher.test(contingency_table_test_LR_flat)

test_that("fisher test is done correctly on LRIs", {
  expect_equal(
    scdiffcom_objects$cond_stat@ora_table$LRI[VALUE == "Apoe:Ldlr"]$OR_UP,
    fisher_LR_up$estimate,
    ignore_attr = TRUE
  )
  expect_equal(
    scdiffcom_objects$cond_stat@ora_table$LRI[VALUE == "Apoe:Ldlr"]$P_VALUE_UP,
    fisher_LR_up$p.value,
    ignore_attr = TRUE
  )
  expect_equal(
    scdiffcom_objects$cond_stat@ora_table$LRI[VALUE == "Csf2:Itgb1"]$OR_DOWN,
    fisher_LR_down$estimate,
    ignore_attr = TRUE
  )
  expect_equal(
    scdiffcom_objects$cond_stat@ora_table$LRI[VALUE == "Csf2:Itgb1"]$P_VALUE_DOWN,
    fisher_LR_down$p.value,
    ignore_attr = TRUE
  )
  expect_equal(
    scdiffcom_objects$cond_stat@ora_table$LRI[VALUE == "Adam10:Axl"]$OR_FLAT,
    fisher_LR_flat$estimate,
    ignore_attr = TRUE
  )
  expect_equal(
    scdiffcom_objects$cond_stat@ora_table$LRI[VALUE == "Adam10:Axl"]$P_VALUE_FLAT,
    fisher_LR_flat$p.value,
    ignore_attr = TRUE
  )
})

# ORA double-check on GO terms

genes_GO_UP <- LRI_mouse$LRI_curated_GO[GO_ID == "GO:0002376"]$LRI

contingency_table_test_GO_up <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_table_detected[LRI %in% genes_GO_UP & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI %in% genes_GO_UP & REGULATION != "UP", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[!(LRI %in% genes_GO_UP) & REGULATION == "UP", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[!(LRI %in% genes_GO_UP) & REGULATION != "UP", .N]
  ),
  2,
  2
)
fisher_GO_up <- fisher.test(contingency_table_test_GO_up)

genes_GO_DOWN <- LRI_mouse$LRI_curated_GO[GO_ID == "GO:0031625"]$LRI
contingency_table_test_GO_down <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_table_detected[LRI %in% genes_GO_DOWN & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI %in% genes_GO_DOWN & REGULATION != "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[!(LRI %in% genes_GO_DOWN) & REGULATION == "DOWN", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[!(LRI %in% genes_GO_DOWN) & REGULATION != "DOWN", .N]
  ),
  2,
  2
)
fisher_GO_down <- fisher.test(contingency_table_test_GO_down)

genes_GO_FLAT <- LRI_mouse$LRI_curated_GO[GO_ID == "GO:0022406"]$LRI
contingency_table_test_GO_flat <- matrix(
  c(
    scdiffcom_objects$cond_stat@cci_table_detected[LRI %in% genes_GO_FLAT & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[LRI %in% genes_GO_FLAT & REGULATION != "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[!(LRI %in% genes_GO_FLAT) & REGULATION == "FLAT", .N],
    scdiffcom_objects$cond_stat@cci_table_detected[!(LRI %in% genes_GO_FLAT) & REGULATION != "FLAT", .N]
  ),
  2,
  2
)
fisher_GO_flat <- fisher.test(contingency_table_test_GO_flat)

test_that(
  "fisher test is done correctly on GO terms",
  {
    expect_equal(
      scdiffcom_objects$cond_stat@ora_table$GO_TERMS[VALUE_BIS == "GO:0002376"]$OR_UP,
      fisher_GO_up$estimate,
      ignore_attr = TRUE
    )
    expect_equal(
      scdiffcom_objects$cond_stat@ora_table$GO_TERMS[VALUE_BIS == "GO:0002376"]$P_VALUE_UP,
      fisher_GO_up$p.value,
      ignore_attr = TRUE
    )
    expect_equal(
      scdiffcom_objects$cond_stat@ora_table$GO_TERMS[VALUE_BIS == "GO:0031625"]$OR_DOWN,
      fisher_GO_down$estimate,
      ignore_attr = TRUE
    )
    expect_equal(
      scdiffcom_objects$cond_stat@ora_table$GO_TERMS[VALUE_BIS == "GO:0031625"]$P_VALUE_DOWN,
      fisher_GO_down$p.value,
      ignore_attr = TRUE
    )
    expect_equal(
      scdiffcom_objects$cond_stat@ora_table$GO_TERMS[VALUE_BIS == "GO:0022406"]$OR_FLAT,
      fisher_GO_flat$estimate,
      ignore_attr = TRUE
    )
    expect_equal(
      scdiffcom_objects$cond_stat@ora_table$GO_TERMS[VALUE_BIS == "GO:0022406"]$P_VALUE_FLAT,
      fisher_GO_flat$p.value,
      ignore_attr = TRUE
    )
  }
)

## Check exported FilterCCI ####

scdiffcom_object_new_filter <- FilterCCI(
  object = scdiffcom_objects$cond_stat,
  new_threshold_quantile_score = 0.25,
  new_threshold_p_value_specificity = 0.01,
  new_threshold_p_value_de = 0.01,
  new_threshold_logfc = 0.01,
  skip_ora = FALSE
)

test_that(
  "cci_table_raw is not modified by FilterCCI",
  {
    expect_identical(
      scdiffcom_objects$cond_stat@cci_table_raw,
      scdiffcom_object_new_filter@cci_table_raw
    )
  }
)

## Check exported RunORA and extra annotations #####

cell_types <-  c(
  "B cell",
  "T cell",
  "endothelial cell of hepatic sinusoid",
  "hepatocyte",
  "myeloid leukocyte"
)
cell_families <- c(
  "leukocyte",
  "leukocyte",
  "endothelial cell",
  "epithelial cell",
  "leukocyte"
)

cell_families_dt <- data.table(
  EMITTER_CELLTYPE = cell_types,
  EMITTER_CELLFAMILY = cell_families
)

scdiffcom_object_new_ora <- RunORA(
  object = scdiffcom_objects$cond_stat,
  extra_annotations = list(
    cell_families_dt
  ),
  overwrite = FALSE
)

test_that(
  "RunORA with extra annotations does not change previous CCI and ORA results",
  {
    expect_identical(
      scdiffcom_objects$cond_stat@cci_table_detected,
      scdiffcom_object_new_ora@cci_table_detected
    )
    expect_identical(
      scdiffcom_objects$cond_stat@ora_table$LRI,
      scdiffcom_object_new_ora@ora_table$LRI
    )
  }
)

## Check accessors ####

retrieved_parameters <- lapply(
  scdiffcom_objects,
  GetParameters
)

test_that(
  "retrieved parameters are the correct ones",
  {
    lapply(
      1:4,
      function(i) {
        expect_identical(
          retrieved_parameters[[i]],
          scdiffcom_objects[[i]]@parameters
        )
      }
    )
  }
)

retrieved_cci_tables <- lapply(
  c(raw = "raw", detected = "detected"),
  function(type) {
    lapply(
      c(simplified = TRUE, non_simplified = FALSE),
      function(simplified) {
        lapply(
          scdiffcom_objects,
          GetTableCCI,
          type = type,
          simplified = simplified
        )
      }
    )
  }
)

test_that(
  "retrieved cci tables are the correct ones",
  {
    expect_identical(
      retrieved_cci_tables$detected$non_simplified$cond_stat,
      scdiffcom_objects$cond_stat@cci_table_detected
    )
  }
)

retrieved_ora_tables <- lapply(
  list("all", c("LRI", "ER_CELLTYPES", "GO_TERMS")),
  function(categories) {
    lapply(
      c(TRUE, FALSE),
      function(simplified) {
        GetTableORA(
          scdiffcom_objects$cond_stat,
          categories = categories,
          simplified = simplified
        )
      }
    )
  }
)

test_that(
  "retrieved ora tables are the correct ones",
  {
    expect_identical(
      retrieved_ora_tables[[2]][[2]]$LRI,
      scdiffcom_objects$cond_stat@ora_table$LRI
    )
  }
)

## Check PlotORA ####

test_that(
  "PlotORA works and returns a ggplot object",
  {
    expect_s3_class(
      PlotORA(
        object = scdiffcom_objects$cond_stat,
        category = "LRI",
        regulation = "DOWN",
        max_terms_show = 50,
        GO_aspect = "biological_process"
      ),
      "ggplot"
    )
    expect_s3_class(
      PlotORA(
        object = scdiffcom_objects$cond_stat,
        category = "ER_CELLTYPES",
        regulation = "DOWN",
        max_terms_show = 20,
        GO_aspect = "biological_process"
      ),
      "ggplot"
    )
  }
)

## Check EraseRawCCI ####

test_that(
  "EraseRawCCI works and returns an empty list for the cci raw table",
  {
    expect_identical(
      EraseRawCCI(scdiffcom_objects$cond_stat)@cci_table_raw,
      list()
    )
    expect_identical(
      EraseRawCCI(scdiffcom_objects$cond_stat)@cci_table_detected,
      scdiffcom_objects$cond_stat@cci_table_detected
    )
  }
)

## Check BuildNetwork ####

types_of_network = c(
  #"condition1_network",
  #"condition2_network",
  #"difference_network",
  #"up_regulated_network",
  #"down_regulated_network",
  "ORA_network"
)

types_of_layout = c(
  "conventional",
  "bipartite"
)

test_that(
  "BuildNetwork works and returns a visNetwork object",
  {
    lapply(
      types_of_network,
      function(net) {
        lapply(
          types_of_layout,
          function(lay) {
            expect_s3_class(
              BuildNetwork(
                object = scdiffcom_objects$cond_stat,
                network_type = net,
                layout_type = lay
              ),
              "visNetwork"
            )
          }
        )
      }
    )
  }
)

## Check BuildNetwork internal functions #####

build_interactive_network(
  object = scdiffcom_objects$cond_stat,
  network_type = types_of_network[[1]],
  layout_type = types_of_layout[[2]],
  class_signature = "scDiffCom",
  subobject_name = NULL,
  abbreviation_table = NULL
)

test_that(
  "No network available for object without statistical analysis",
  {
    expect_error(
      build_interactive_network(
        object = scdiffcom_objects$cond_nostat,
        network_type = types_of_network[[1]],
        layout_type = types_of_layout[[2]],
        class_signature = "scDiffCom",
        subobject_name = NULL,
        abbreviation_table = NULL
      )
    )
  }
)
