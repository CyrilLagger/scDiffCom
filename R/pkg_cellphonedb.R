#' Create text files to be used by CellPhoneDB from a Seurat object
#'
#' @param seurat_obj a Seurat object
#' @param assay assay to pull data from
#' @param slot slot to pull data from
#' @param log_scale logical
#' @param seurat_cell_type_id xx
#' @param min_cells xx
#' @param condition_id xx
#' @param input_dir xx
#'
#' @return Path where the data and metadata are written.
#' @export
#'
#' @examples
#'
create_cpdb_input <- function(seurat_obj,
                              assay = "RNA",
                              slot = "data",
                              log_scale = TRUE,
                              seurat_cell_type_id,
                              min_cells = 10,
                              condition_id = NULL,
                              input_dir = getwd()
) {
  data <- prepare_seurat_data(seurat_obj = seurat_obj,
                              assay = assay,
                              slot = slot,
                              log_scale = log_scale,
                              convert_to_human = TRUE,
                              return_type = "data.frame")
  metadata <- prepare_seurat_metadata(seurat_obj = seurat_obj,
                                      seurat_cell_type_id = seurat_cell_type_id,
                                      condition_id = condition_id)
  cell_type_filt <- filter_cell_types(metadata = metadata,
                                      min_cells = min_cells)
  metadata <- metadata[metadata$cell_type %in% cell_type_filt, ]
  data <- data[, colnames(data) %in% metadata$cell_id]
  data <- tibble::rownames_to_column(data, var = "Gene")
  if(is.null(condition_id)) {
    colnames(metadata) <- c("Cell", "cell_type")
    message(paste0("Writing CellPhoneDB input data to ", input_dir, "/cpdb_data_noCond.txt"))
    utils::write.table(data,
                       file = paste0(input_dir, "/cpdb_data_noCond.txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')
    message(paste0("Writing CellPhoneDB input metadata to ", input_dir, "/cpdb_metadata_noCond.txt"))
    utils::write.table(metadata,
                       file = paste0(input_dir, "/cpdb_metadata_noCond.txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')
    return(list(data_path = paste0(input_dir, "/cpdb_data_noCond.txt"),
                metadata_path = paste0(input_dir, "/cpdb_metadata_noCond.txt") ))
  } else {
    conds <- unique(metadata$condition)
    if(length(conds) != 2) stop("Wrong number of groups in cell-type conditions (expected 2).")
    meta1 <- metadata[metadata$condition == conds[[1]], ]
    data1 <- data[, colnames(data) %in% c("Gene", meta1$cell_id)]
    meta1$condition <- NULL
    colnames(meta1) <- c("Cell", "cell_type")
    meta2 <- metadata[metadata$condition == conds[[2]], ]
    data2 <- data[, colnames(data) %in% c("Gene", meta2$cell_id)]
    meta2$condition <- NULL
    colnames(meta2) <- c("Cell", "cell_type")
    message(paste0("Writing CellPhoneDB input data to ", input_dir, "/cpdb_data_" , conds[[1]], ".txt"))
    utils::write.table(data1,
                       file = paste0(input_dir, "/cpdb_data_" , conds[[1]], ".txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')
    message(paste0("Writing CellPhoneDB input data to ", input_dir, "/cpdb_data_" , conds[[2]], ".txt"))
    utils::write.table(data2,
                       file = paste0(input_dir, "/cpdb_data_" , conds[[2]], ".txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')
    message(paste0("Writing CellPhoneDB input metadata to ", input_dir, "/cpdb_metadata_", conds[[1]], ".txt"))
    utils::write.table(meta1,
                       file = paste0(input_dir, "/cpdb_metadata_", conds[[1]], ".txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')
    message(paste0("Writing CellPhoneDB input metadata to ", input_dir, "/cpdb_metadata_", conds[[2]], ".txt"))
    utils::write.table(metadata,
                       file = paste0(input_dir, "/cpdb_metadata_", conds[[2]], ".txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')

    return(list(data_path1 = paste0(input_dir, "/cpdb_data_" , conds[[1]], ".txt"),
                metadata_path1 = paste0(input_dir, "/cpdb_metadata_", conds[[1]], ".txt"),
                data_path2 = paste0(input_dir, "/cpdb_data_" , conds[[2]], ".txt"),
                metadata_path2 = paste0(input_dir, "/cpdb_metadata_", conds[[2]], ".txt")))
  }
}


#' Run CellPhoneDB from text files.
#'
#' @param data_path a text file
#' @param metadata_path a text file
#' @param method character
#' @param project_name character
#' @param iterations integer
#' @param threshold numeric
#' @param result_precision numeric
#' @param counts_data character, type of gene identifier
#' @param output_path character
#' @param output_format character
#' @param means_result_name character
#' @param significant_mean_result_name character
#' @param deconvoluted_result_name character
#' @param verbose logical
#' @param pvalues_result_name character
#' @param debug_seed numeric
#' @param threads numeric
#' @param subsampling logical
#' @param subsampling_log logical
#' @param subsampling_num_pc numeric
#' @param subsampling_num_cells numeric
#'
#' @return Run CellPhoneDB and return results in output directory
#' @export
#'
#' @examples
#'
run_cpdb_from_files <- function(data_path,
                                metadata_path,
                                method = 'statistical_analysis',
                                project_name = NULL,
                                iterations = NULL,
                                threshold = NULL,
                                result_precision = NULL,
                                counts_data = NULL,
                                output_path = NULL,
                                output_format = NULL,
                                means_result_name = NULL,
                                significant_mean_result_name = NULL,
                                deconvoluted_result_name = NULL,
                                verbose = TRUE,
                                pvalues_result_name = NULL,
                                debug_seed = NULL,
                                threads = NULL,
                                subsampling = FALSE,
                                subsampling_log = FALSE,
                                subsampling_num_pc = NULL,
                                subsampling_num_cells = NULL
) {
  command_cpdb <- "cellphonedb method"
  command_cpdb <- paste(command_cpdb, method, metadata_path, data_path, sep = " ")
  if(method == "statistical_analysis") {
    if(subsampling) {
      command_cpdb <- paste0(command_cpdb, " --subsampling --subsampling-log")
      if(subsampling_log) {
        command_cpdb <- paste0(command_cpdb, " true")
      } else {
        command_cpdb <- paste0(command_cpdb, " false")
      }
      if (!is.null(subsampling_num_pc)) command_cpdb <- paste0(command_cpdb, " --subsampling-num-pc=", subsampling_num_pc)
      if (!is.null(subsampling_num_cells)) command_cpdb <- paste0(command_cpdb, " --subsampling-num-cells=", subsampling_num_cells)
    }
    if (!is.null(pvalues_result_name)) command_cpdb <- paste0(command_cpdb, " --pvalues-result-name=", pvalues_result_name)
    if (!is.null(debug_seed)) command_cpdb <- paste0(command_cpdb, " --debug-seed=", debug_seed)
    if (!is.null(threads)) command_cpdb <- paste0(command_cpdb, " --threads=", threads)
  }
  if (!is.null(project_name)) command_cpdb <- paste0(command_cpdb, " --project-name=", project_name)
  if (!is.null(iterations)) command_cpdb <- paste0(command_cpdb, " --iterations=", iterations)
  if (!is.null(result_precision)) command_cpdb <- paste0(command_cpdb, " --result-precision=", result_precision)
  if (!is.null(counts_data)) command_cpdb <- paste0(command_cpdb, " --counts-data=", counts_data)
  if (!is.null(output_path)) command_cpdb <- paste0(command_cpdb, " --output-path=", output_path)
  if (!is.null(output_format)) command_cpdb <- paste0(command_cpdb, " --output-format=", output_format)
  if (!is.null(means_result_name)) command_cpdb <- paste0(command_cpdb, " --means-result-name=", means_result_name)
  if (!is.null(significant_mean_result_name)) command_cpdb <- paste0(command_cpdb, " --significant-mean-result-name=", significant_mean_result_name)
  if (!is.null(deconvoluted_result_name)) command_cpdb <- paste0(command_cpdb, " --deconvoluted-result-name=", deconvoluted_result_name)
  if(verbose) {
    command_cpdb <- paste0(command_cpdb, " --verbose")
  } else {
    command_cpdb <- paste0(command_cpdb, " --quiet")
  }
  system(command = command_cpdb)
}

#' Run CellPhoneDB from a Seurat object
#'
#' @param seurat_obj x
#' @param assay x
#' @param slot x
#' @param log_scale x
#' @param seurat_cell_type_id x
#' @param condition_id x
#' @param min_cells x
#' @param input_dir x
#' @param method x
#' @param iterations x
#' @param threshold x
#' @param result_precision x
#' @param counts_data x
#' @param output_format x
#' @param means_result_name x
#' @param significant_mean_result_name x
#' @param deconvoluted_result_name x
#' @param verbose x
#' @param pvalues_result_name x
#' @param debug_seed x
#' @param threads x
#' @param subsampling x
#' @param subsampling_log x
#' @param subsampling_num_pc x
#' @param subsampling_num_cells x
#'
#' @return Return the results of CellPhoneDB
#' @export
#'
#' @examples
#'
run_cpdb_from_seurat <- function(seurat_obj,
                                 assay = "RNA",
                                 slot = "data",
                                 log_scale = TRUE,
                                 seurat_cell_type_id,
                                 min_cells = 10,
                                 condition_id = NULL,
                                 input_dir = getwd(),
                                 method = 'statistical_analysis',
                                 iterations = NULL,
                                 threshold = NULL,
                                 result_precision = NULL,
                                 counts_data = NULL,
                                 output_format = NULL,
                                 means_result_name = NULL,
                                 significant_mean_result_name = NULL,
                                 deconvoluted_result_name = NULL,
                                 verbose = TRUE,
                                 pvalues_result_name = NULL,
                                 debug_seed = NULL,
                                 threads = NULL,
                                 subsampling = FALSE,
                                 subsampling_log = FALSE,
                                 subsampling_num_pc = NULL,
                                 subsampling_num_cells = NULL
) {
  paths <- create_cpdb_input(seurat_obj = seurat_obj,
                             assay = assay,
                             slot = slot,
                             log_scale = log_scale,
                             seurat_cell_type_id = seurat_cell_type_id,
                             condition_id = condition_id,
                             min_cells = min_cells,
                             input_dir = input_dir)
  if(is.null(condition_id)) {
    run_cpdb_from_files(data_path = paths$data_path,
                        metadata_path = paths$metadata_path,
                        method = method,
                        project_name = "cpdb_results_noCond",
                        iterations = iterations,
                        threshold = threshold,
                        result_precision = result_precision,
                        counts_data = counts_data,
                        output_path = input_dir,
                        output_format = output_format,
                        means_result_name = means_result_name,
                        significant_mean_result_name = significant_mean_result_name,
                        deconvoluted_result_name = deconvoluted_result_name,
                        verbose = verbose,
                        pvalues_result_name = pvalues_result_name,
                        debug_seed = debug_seed,
                        threads = threads,
                        subsampling = subsampling,
                        subsampling_log = subsampling_log,
                        subsampling_num_pc = subsampling_num_pc,
                        subsampling_num_cells = subsampling_num_cells
    )
  } else {
    run_cpdb_from_files(data_path = paths$data_path1,
                        metadata_path = paths$metadata_path1,
                        method = method,
                        project_name = "cpdb_results_Cond1",
                        iterations = iterations,
                        threshold = threshold,
                        result_precision = result_precision,
                        counts_data = counts_data,
                        output_path = input_dir,
                        output_format = output_format,
                        means_result_name = means_result_name,
                        significant_mean_result_name = significant_mean_result_name,
                        deconvoluted_result_name = deconvoluted_result_name,
                        verbose = verbose,
                        pvalues_result_name = pvalues_result_name,
                        debug_seed = debug_seed,
                        threads = threads,
                        subsampling = subsampling,
                        subsampling_log = subsampling_log,
                        subsampling_num_pc = subsampling_num_pc,
                        subsampling_num_cells = subsampling_num_cells
    )
    run_cpdb_from_files(data_path = paths$data_path2,
                        metadata_path = paths$metadata_path2,
                        method = method,
                        project_name = "cpdb_results_Cond2",
                        iterations = iterations,
                        threshold = threshold,
                        result_precision = result_precision,
                        counts_data = counts_data,
                        output_path = input_dir,
                        output_format = output_format,
                        means_result_name = means_result_name,
                        significant_mean_result_name = significant_mean_result_name,
                        deconvoluted_result_name = deconvoluted_result_name,
                        verbose = verbose,
                        pvalues_result_name = pvalues_result_name,
                        debug_seed = debug_seed,
                        threads = threads,
                        subsampling = subsampling,
                        subsampling_log = subsampling_log,
                        subsampling_num_pc = subsampling_num_pc,
                        subsampling_num_cells = subsampling_num_cells
    )
  }
}
