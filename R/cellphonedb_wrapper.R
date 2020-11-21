#' Run CellPhoneDB from a Seurat object
#'
#' @param seurat_object x
#' @param assay x
#' @param slot x
#' @param log_scale x
#' @param celltype_column_id x
#' @param input_species x
#' @param condition_column_id x
#' @param min_cells x
#' @param input_dir x
#' @param create_plots x
#' @param method x
#' @param iterations x
#' @param threshold x
#' @param result_precision x
#' @param counts_data x
#' @param output_format x
#' @param verbose x
#' @param debug_seed x
#' @param threads x
#' @param subsampling x
#' @param subsampling_log x
#' @param subsampling_num_pc x
#' @param subsampling_num_cells x
#' @param return_full_dt x
#'
#' @return Return the results of CellPhoneDB
#' @export
run_cpdb_from_seurat <- function(
                                 seurat_object,
                                 assay = "RNA",
                                 slot = "data",
                                 log_scale = FALSE,
                                 celltype_column_id,
                                 input_species = "mouse",
                                 min_cells = 5,
                                 condition_column_id = NULL,
                                 input_dir = getwd(),
                                 create_plots = FALSE,
                                 method = "statistical_analysis",
                                 iterations = NULL,
                                 threshold = NULL,
                                 result_precision = NULL,
                                 counts_data = NULL,
                                 output_format = NULL,
                                 verbose = TRUE,
                                 debug_seed = NULL,
                                 threads = NULL,
                                 subsampling = FALSE,
                                 subsampling_log = FALSE,
                                 subsampling_num_pc = NULL,
                                 subsampling_num_cells = NULL,
                                 return_full_dt = TRUE) {
  message("Create file directory if not already existing.")
  if (!dir.exists(input_dir)) {
    dir.create(input_dir)
  }
  message("Create input files from Seurat to be used by CellPhoneDB.")
  paths <- create_cpdb_input(
    seurat_object = seurat_object,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    celltype_column_id = celltype_column_id,
    input_species = input_species,
    condition_column_id = condition_column_id,
    min_cells = min_cells,
    input_dir = input_dir
  )
  if (is.null(condition_column_id)) {
    n_run <- 1
    project_name <- "cpdb_results_noCond"
    means_result_name <- NULL
    significant_means_result_name <- NULL
    econvoluted_result_name <- NULL
    pvalues_result_name <- NULL
    message("Run CellphoneDB once (no condition).")
  } else {
    n_run <- 2
    project_name <- paste0("cpdb_results_", paths$conds)
    means_result_name <- paste0("means-", paths$conds, ".txt")
    significant_means_result_name <- paste0("significant-means-", paths$conds, ".txt")
    deconvoluted_result_name <- paste0("deconvoluted-", paths$conds, ".txt")
    pvalues_result_name <- paste0("pvalues-", paths$conds, ".txt")
    message("Run CellPhoneDB on the two conditions.")
  }
  for (i in 1:n_run) {
    run_cpdb_from_files(
      data_path = paths$data_path[[i]],
      metadata_path = paths$md_path[[i]],
      method = method,
      project_name = project_name[[i]],
      iterations = iterations,
      threshold = threshold,
      result_precision = result_precision,
      counts_data = counts_data,
      output_path = input_dir,
      output_format = output_format,
      means_result_name = means_result_name[[i]],
      significant_means_result_name = significant_means_result_name[[i]],
      deconvoluted_result_name = deconvoluted_result_name[[i]],
      verbose = verbose,
      pvalues_result_name = pvalues_result_name[[i]],
      debug_seed = debug_seed,
      threads = threads,
      subsampling = subsampling,
      subsampling_log = subsampling_log,
      subsampling_num_pc = subsampling_num_pc,
      subsampling_num_cells = subsampling_num_cells
    )
  }
  if (return_full_dt) {
    if (method != "statistical_analysis") {
      message("Not possible to create the full data.table for the selected method.")
    } else {
      message("Create and write full data.table.")
      full_dt <- create_cpdb_cci(
        input_dir = input_dir,
        conds = paths$conds
      )
      if (is.null(condition_column_id)) {
        output_dir_full <- paste0(input_dir, "/cpdb_full_table_noCond.txt")
      } else {
        output_dir_full <- paste0(input_dir, "/cpdb_full_table_withCond.txt")
      }
      utils::write.table(
        full_dt,
        file = output_dir_full,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = "\t"
      )
    }
  }
  if (input_species == "mouse") {
    message(paste0("Write human-mouse orthologs in ", input_dir, "/cpdb_human_mouse_orthologs.txt"))
    utils::write.table(
      paths$gene_mapping,
      file = paste0(input_dir, "/cpdb_human_mouse_orthologs.txt"),
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }
}
