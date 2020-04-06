#' Create text files to be used by CellPhoneDB from a Seurat object
#'
#' @param seurat_obj a Seurat object
#' @param assay assay to pull data from
#' @param slot slot to pull data from
#' @param log_scale logical
#' @param seurat_cell_type_id xx
#' @param condition_id xx
#' @param input_dir xx
#'
#' @return none
#' @export
#'
#' @examples
#'
create_cpdb_input <- function(seurat_obj,
                                  assay = "RNA",
                                  slot = "data",
                                  log_scale = TRUE,
                                  seurat_cell_type_id,
                                  condition_id = NULL,
                                  input_dir = getwd()
) {
  data <- prepare_seurat_data(seurat_obj = seurat_obj,
                              assay = assay,
                              slot = slot,
                              log_scale = log_scale,
                              convert_to_human = TRUE,
                              return_type = "data.frame")
  if(is.null(condition_id)) {
    metadata <- prepare_seurat_metadata(seurat_obj = seurat_obj,
                                        seurat_cell_type_id = seurat_cell_type_id,
                                        condition_id = NULL)
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
  } else {
    stop("Applying CellPhoneDB to multiple conditions is not yet implemented. Stay tuned!")
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
#' @param input_dir x
#' @param method x
#' @param project_name x
#' @param iterations x
#' @param threshold x
#' @param result_precision x
#' @param counts_data x
#' @param output_path x
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
                                 condition_id = NULL,
                                 input_dir = getwd(),
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
  create_cpdb_input(seurat_obj = seurat_obj,
                    assay = assay,
                    slot = slot,
                    log_scale = log_scale,
                    seurat_cell_type_id = seurat_cell_type_id,
                    condition_id = condition_id,
                    input_dir = input_dir)
  if(is.null(condition_id)) {
    data_path <- paste0(input_dir, "/cpdb_data_noCond.txt")
    metadata_path <- paste0(input_dir, "/cpdb_metadata_noCond.txt")
    run_cpdb_from_files(data_path = data_path,
                        metadata_path = metadata_path,
                        method = method,
                        project_name = project_name,
                        iterations = iterations,
                        threshold = threshold,
                        result_precision = result_precision,
                        counts_data = counts_data,
                        output_path = output_path,
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
    stop("Applying CellPhoneDB to multiple conditions is not yet implemented. Stay tuned!")
  }
}
