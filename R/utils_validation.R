#' Check the parameters
#'
#' @param params A list of parameters.
#'
#' @noRd
validate_parameters <- function(
  params,
  from_inputs
) {
  params_names_base <- c(
    "object_name",
    "LRdb_species",
    "seurat_celltype_id",
    "seurat_condition_id",
    "cond1_name",
    "cond2_name",
    "seurat_assay",
    "seurat_slot",
    "log_scale",
    "threshold_min_cells",
    "threshold_pct",
    "permutation_analysis",
    "iterations",
    "threshold_quantile_score",
    "threshold_p_value_specificity",
    "threshold_p_value_de",
    "threshold_logfc",
    "return_distributions",
    "seed",
    "verbose"
  )
  params_names_additional <- c(
    "conditional_analysis",
    "max_nL",
    "max_nR"
  )
  params_names_all <- c(params_names_base, params_names_additional)
  if (from_inputs) {
    if (!identical(sort(names(params)), sort(params_names_base))) {
      stop("Parameters do not match.")
    }
    params <- params[params_names_base]
  } else {
    if (!identical(sort(names(params)), sort(params_names_all))) {
      stop("Parameters do not match.")
    }
    params <- params[params_names_all]
  }
  res <- NULL
  if (!is.character(params$object_name) | length(params$object_name) != 1) {
    res <- c(res, "`object_name` must be a character vector of length 1")
  }
  if (!(params$LRdb_species %in% c("mouse", "human"))) {
    res <- c(res, "`LRdb_species` must be either 'human' or 'mouse'")
  }
  if (!is.character(params$seurat_celltype_id) | length(params$seurat_celltype_id) != 1) {
    res <- c(res, "`seurat_celltype_id` must be a character vector of length 1")
  }
  if (!is.null(params$seurat_condition_id)) {
    if (!is.character(params$seurat_condition_id) | length(params$seurat_condition_id) != 1) {
      res <- c(res, "`seurat_condition_id` must be NULL or a character vector of length 1")
    }
  }
  if (!is.null(params$cond1_name)) {
    if (!is.character(params$cond1_name) | length(params$cond1_name) != 1) {
      res <- c(res, "`cond1_name` must be NULL or a character vector of length 1")
    }
  }
  if (!is.null(params$cond2_name)) {
    if (!is.character(params$cond2_name) | length(params$cond2_name) != 1) {
      res <- c(res, "`cond2_name` must be NULL or a character vector of length 1")
    }
  }
  if (!is.character(params$seurat_assay) | length(params$seurat_assay) != 1) {
    res <- c(res, "`seurat_assay` must be NULL or a character vector of length 1")
  }
  if (!(params$seurat_slot %in% c("counts", "data"))) {
    res <- c(res, "`seurat_slot` must be either `data` or `counts`")
  }
  if(!is.logical(params$log_scale) | length(params$log_scale) != 1) {
    res <- c(res, "`log_scale` must be a logical vector of length 1")
  }
  if (!is.numeric(params$threshold_min_cells) | length(params$threshold_min_cells) > 1) {
    res <- c(res, "`threshold_min_cells` must be a numeric vector of length 1")
  } else if (params$threshold_min_cells < 0) {
    res <- c(res, "`threshold_min_cells` must be a non-negative numeric")
  }
  if (!is.numeric(params$threshold_pct) | length(params$threshold_pct) != 1) {
    res <- c(res, "`threshold_pct` must be a numeric vector of length 1")
  } else if(params$threshold_pct < 0 | params$threshold_pct >= 1) {
    res <- c(res, "`threshold_pct` must be a numeric in [0,1[")
  }
  if(!is.logical(params$permutation_analysis) | length(params$permutation_analysis) != 1) {
    res <- c(res, "`permutation_analysis` must be a logical vector of length 1")
  }
  if (!is.numeric(params$iterations) | length(params$iterations) > 1) {
    res <- c(res, "`iterations` must be a numeric vector of length 1")
  } else if (params$iterations <= 0) {
    res <- c(res, "`iterations` must be a positive integer")
  } else  if (params$iterations > 20000) {
    warning("performing more than 20'000 `iterations` will require a lot of memory")
  }
  if (!is.numeric(params$threshold_quantile_score) | length(params$threshold_quantile_score) != 1) {
    res <- c(res, "`threshold_quantile_score` must be a numeric vector of length 1")
  } else if (params$threshold_quantile_score < 0 |params$ threshold_quantile_score >= 1) {
    res <- c(res, "`threshold_quantile_score` must be a numeric in [0,1[")
  }
  if (!is.numeric(params$threshold_p_value_specificity) | length(params$threshold_p_value_specificity) != 1) {
    res <- c(res, "`threshold_p_value_specificity` must be a numeric vector of length 1")
  } else if (params$threshold_p_value_specificity <= 0 | params$threshold_p_value_specificity > 1) {
    res <- c(res, "`threshold_p_value_specificity` must be a numeric in ]0,1]")
  }
  if (!is.numeric(params$threshold_p_value_de) | length(params$threshold_p_value_de) != 1) {
    res <- c(res, "`threshold_p_value_de` must be a numeric vector of length 1")
  } else if(params$threshold_p_value_de <= 0 | params$threshold_p_value_de > 1) {
    res <- c(res, "`threshold_p_value_de` must be a numeric in ]0,1]")
  }
  if (!is.numeric(params$threshold_logfc) | length(params$threshold_logfc) != 1) {
    res <- c(res, "`threshold_logfc` must be a numeric vector of length 1")
  } else if(params$threshold_logfc <= 0) {
    res <- c(res, "`threshold_logfc` must be a positive numeric")
  }
  if(!is.logical(params$return_distributions) | length(params$return_distributions) != 1) {
    res <- c(res, "`return_distributions` must be a logical vector of length 1")
  }
  if(!is.logical(params$verbose) | length(params$verbose) != 1) {
    res <- c(res, "`verbose` must be a logical vector of length 1")
  }
  if(!from_inputs) {
    if(!is.logical(params$conditional_analysis) | length(params$conditional_analysis) != 1) {
      res <- c(res, "`conditional_analysis` must be a logical vector of length 1")
    }
    if (!is.numeric(params$max_nL) | length(params$max_nL) != 1) {
      res <- c(res, "`max_nL` must be a numeric vector of length 1")
    }
    if (!is.numeric(params$max_nR) | length(params$max_nR) != 1) {
      res <- c(res, "`max_nR` must be a numeric vector of length 1")
    }
  }
  list(
    params = params,
    check =res
  )
}

#' Check the validity of \code{parameters}
#'
#' @param parameters A list of parameters to check the validity of.
#'
#' @noRd
validate_slot_parameters <- function(parameters) {
  res <- validate_parameters(
    params = parameters,
    from_inputs = FALSE
  )$check
  if(is.null(res)){
    NULL
  } else {
    "@parameters is not formatted the correct way"
  }
}

#' Check the validity of \code{cci_raw}
#'
#' @param parameters A list of parameters.
#' @param cci_raw A data.table to check the validity of.
#'
#' @noRd
validate_slot_cci_raw <- function(
  parameters,
  cci_raw
) {
  NULL
}

#' Check the validity of \code{cci_detected}
#'
#' @param parameters A list of parameters.
#' @param cci_detected A data.table to check the validity of.
#'
#' @noRd
validate_slot_cci_detected <- function(
  parameters,
  cci_detected
) {
  NULL
}

#' Check the validity of \code{ora_default}
#'
#' @param parameters A list of parameters.
#' @param ora_default A list of data.tables to check the validity of.
#'
#' @noRd
validate_slot_ora_default <- function(
  parameters,
  ora_default
) {
  NULL
}

#' Check the validity of \code{ora_stringent}
#'
#' @param parameters A list of parameters.
#' @param ora_default A list of data.tables to check the validity of.
#'
#' @noRd
validate_slot_ora_stringent <- function(
  parameters,
  ora_default
) {
  NULL
}

#' Check the validity of \code{distributions}
#'
#' @param parameters A list of parameters.
#' @param distributions A list of matrices to check the validity of.
#'
#' @noRd
validate_slot_distributions <- function(
  parameters,
  distributions
) {
  NULL
}

validate_slot_is_combined <- function()
{
  NULL
}
