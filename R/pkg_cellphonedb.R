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
  if(is.null(condition_id)) {
    data <- tibble::rownames_to_column(data, var = "Gene")
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
    data1 <- data[, colnames(data) %in% meta1$cell_id]
    data1 <- tibble::rownames_to_column(data1, var = "Gene")
    meta1$condition <- NULL
    colnames(meta1) <- c("Cell", "cell_type")
    meta2 <- metadata[metadata$condition == conds[[2]], ]
    data2 <- data[, colnames(data) %in% meta2$cell_id]
    data2 <- tibble::rownames_to_column(data2, var = "Gene")
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
    utils::write.table(meta2,
                       file = paste0(input_dir, "/cpdb_metadata_", conds[[2]], ".txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')

    return(list(cond1 = conds[[1]],
                cond2 = conds[[2]],
                data_path1 = paste0(input_dir, "/cpdb_data_" , conds[[1]], ".txt"),
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
#' @param significant_means_result_name character
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
                                significant_means_result_name = NULL,
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
  if (!is.null(significant_means_result_name)) command_cpdb <- paste0(command_cpdb, " --significant-means-result-name=", significant_means_result_name)
  if (!is.null(deconvoluted_result_name)) command_cpdb <- paste0(command_cpdb, " --deconvoluted-result-name=", deconvoluted_result_name)
  if(verbose) {
    command_cpdb <- paste0(command_cpdb, " --verbose")
  } else {
    command_cpdb <- paste0(command_cpdb, " --quiet")
  }
  system(command = command_cpdb)
}

#' Title
#'
#' @param means_path x
#' @param pvalues_path x
#' @param output_path x
#' @param output_name x
#' @param rows x
#' @param columns x
#' @param verbose x
#'
#' @return x
#' @export
#'
#' @examples
run_cpdb_dot_plot <- function(means_path = NULL,
                              pvalues_path = NULL,
                              output_path = NULL,
                              output_name = NULL,
                              rows = NULL,
                              columns = NULL,
                              verbose = TRUE

) {
  command_dp <- "cellphonedb plot dot_plot"

  if (!is.null(means_path)) command_dp <- paste0(command_dp, " --means-path=", means_path)
  if (!is.null(pvalues_path)) command_dp <- paste0(command_dp, " --pvalues-path=", pvalues_path)
  if (!is.null(output_path)) command_dp <- paste0(command_dp, " --output-path=", output_path)
  if (!is.null(output_name)) command_dp <- paste0(command_dp, " --output-name=", output_name)
  if (!is.null(rows)) command_dp <- paste0(command_dp, " --rows=", rows)
  if (!is.null(columns)) command_dp <- paste0(command_dp, " --columns=", columns)
  if(verbose) {
    command_dp <- paste0(command_dp, " --verbose")
  } else {
    command_dp <- paste0(command_dp, " --quiet")
  }
  system(command = command_dp)
}

#' Title
#'
#' @param metadata_path x
#' @param pvalues_path x
#' @param output_path x
#' @param count_name x
#' @param log_name x
#' @param count_network_name x
#' @param interaction_count_name x
#' @param verbose x
#'
#' @return x
#' @export
#'
#' @examples
run_cpdb_heatmap <- function(metadata_path,
                             pvalues_path = NULL,
                             output_path = NULL,
                             count_name = NULL,
                             log_name = NULL,
                             count_network_name = NULL,
                             interaction_count_name = NULL,
                             verbose = TRUE

) {
  command_hm <- "cellphonedb plot heatmap_plot"
  command_hm <- paste(command_hm, metadata_path, sep = " ")

  if (!is.null(pvalues_path)) command_hm <- paste0(command_hm, " --pvalues-path=", pvalues_path)
  if (!is.null(output_path)) command_hm <- paste0(command_hm, " --output-path=", output_path)
  if (!is.null(count_name)) command_hm <- paste0(command_hm, " --count-name=", count_name)
  if (!is.null(log_name)) command_hm <- paste0(command_hm, " --log-name=", log_name)
  if (!is.null(count_network_name)) command_hm <- paste0(command_hm, " --count-network-name=", count_network_name)
  if (!is.null(interaction_count_name)) command_hm <- paste0(command_hm, " --interaction-count-name=", interaction_count_name)
  if(verbose) {
    command_hm <- paste0(command_hm, " --verbose")
  } else {
    command_hm <- paste0(command_hm, " --quiet")
  }
  system(command = command_hm)

}

#' Title
#'
#' @param input_dir x
#' @param condition_id x
#' @param cond1 x
#' @param cond2 x
#'
#' @return x
#' @export
#'
#' @examples
create_cpdp_cci <- function(
  input_dir,
  condition_id = NULL,
  cond1 = NULL,
  cond2 = NULL
) {
  if(is.null(condition_id)) {
    cpdb_means <- utils::read.table(file = paste0(input_dir, '/cpdb_results_noCond/means.txt'),
                                    header = TRUE,
                                    sep = "\t")
    cpdb_pvalues <- utils::read.table(file = paste0(input_dir, '/cpdb_results_noCond/pvalues.txt'),
                                      header = TRUE,
                                      sep = "\t")
    if(!identical(cpdb_means[, 1:11], cpdb_pvalues[, 1:11]) |
       !identical(colnames(cpdb_means), colnames(cpdb_pvalues))) {
      stop("Non identical columns or interactions in means.txt and pvalues.txt.")
    }
    long_means <- data.table::melt(setDT(cpdb_means), id.vars = colnames(cpdb_means)[1:11], variable.name = "cell-type_pair", value.name = "score")
    long_pvalues <- data.table::melt(setDT(cpdb_pvalues), id.vars = colnames(cpdb_pvalues)[1:11], variable.name = "cell-type_pair", value.name = "pvalue")
    cpdb_comb <- data.table::merge.data.table(long_means, long_pvalues)
    return(cpdb_comb)
  } else {
    cpdb_means_cond1 <- utils::read.table(file = paste0(input_dir, "/cpdb_results_", cond1, "/means-", cond1, ".txt"),
                                    header = TRUE,
                                    sep = "\t")
    cpdb_pvalues_cond1 <- utils::read.table(file = paste0(input_dir, "/cpdb_results_", cond1, "/pvalues-", cond1, ".txt"),
                                      header = TRUE,
                                      sep = "\t")
    cpdb_means_cond2 <- utils::read.table(file = paste0(input_dir, "/cpdb_results_", cond2, "/means-", cond2, ".txt"),
                                          header = TRUE,
                                          sep = "\t")
    cpdb_pvalues_cond2 <- utils::read.table(file = paste0(input_dir, "/cpdb_results_", cond2, "/pvalues-", cond2, ".txt"),
                                            header = TRUE,
                                            sep = "\t")
    if(!identical(cpdb_means_cond1[, 1:11], cpdb_pvalues_cond1[, 1:11]) |
       !identical(colnames(cpdb_means_cond1), colnames(cpdb_pvalues_cond1)) |
       !identical(cpdb_means_cond2[, 1:11], cpdb_pvalues_cond2[, 1:11]) |
       !identical(colnames(cpdb_means_cond2), colnames(cpdb_pvalues_cond2))) {
      stop("Non identical columns or interactions in means and pvalues files.")
    }
    long_means_cond1 <- data.table::melt(setDT(cpdb_means_cond1), id.vars = colnames(cpdb_means_cond1)[1:11], variable.name = "cell-type_pair", value.name = "score")
    long_pvalues_cond1 <- data.table::melt(setDT(cpdb_pvalues_cond1), id.vars = colnames(cpdb_pvalues_cond1)[1:11], variable.name = "cell-type_pair", value.name = "pvalue")
    cpdb_comb_cond1 <- data.table::merge.data.table(long_means_cond1, long_pvalues_cond1)
    long_means_cond2 <- data.table::melt(setDT(cpdb_means_cond2), id.vars = colnames(cpdb_means_cond2)[1:11], variable.name = "cell-type_pair", value.name = "score")
    long_pvalues_cond2 <- data.table::melt(setDT(cpdb_pvalues_cond2), id.vars = colnames(cpdb_pvalues_cond2)[1:11], variable.name = "cell-type_pair", value.name = "pvalue")
    cpdb_comb_cond2 <- data.table::merge.data.table(long_means_cond2, long_pvalues_cond2)

    cpdb_comb <- merge.data.table(cpdb_comb_cond1, cpdb_comb_cond2, suffixes = c(paste0("_", cond1), paste0("_", cond2)), all = TRUE)
    return(cpdb_comb)
  }
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
#' @param verbose x
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
                                 verbose = TRUE,
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
                        means_result_name = NULL,
                        significant_means_result_name = NULL,
                        deconvoluted_result_name = NULL,
                        verbose = verbose,
                        pvalues_result_name = NULL,
                        debug_seed = debug_seed,
                        threads = threads,
                        subsampling = subsampling,
                        subsampling_log = subsampling_log,
                        subsampling_num_pc = subsampling_num_pc,
                        subsampling_num_cells = subsampling_num_cells
    )
    run_cpdb_dot_plot(means_path = paste0(input_dir, "/cpdb_results_noCond/means.txt"),
                      pvalues_path = paste0(input_dir, "/cpdb_results_noCond/pvalues.txt"),
                      output_path = paste0(input_dir, "/cpdb_results_noCond"),
                      output_name = "dot_plot.pdf",
                      rows = NULL,
                      columns = NULL,
                      verbose = verbose)
    run_cpdb_heatmap(metadata_path = paths$metadata_path,
                     pvalues_path = paste0(input_dir, "/cpdb_results_noCond/pvalues.txt"),
                     output_path = paste0(input_dir, "/cpdb_results_noCond"),
                     count_name = "heatmap_count.pdf",
                     log_name = "heatmap_log_count.pdf",
                     count_network_name = "network.txt",
                     interaction_count_name = "interaction_count.txt",
                     verbose = verbose)
  } else {
    run_cpdb_from_files(data_path = paths$data_path1,
                        metadata_path = paths$metadata_path1,
                        method = method,
                        project_name = paste0("cpdb_results_", paths$cond1),
                        iterations = iterations,
                        threshold = threshold,
                        result_precision = result_precision,
                        counts_data = counts_data,
                        output_path = input_dir,
                        output_format = output_format,
                        means_result_name = paste0("means-", paths$cond1, ".txt"),
                        significant_means_result_name = paste0("significant-means-", paths$cond1, ".txt"),
                        deconvoluted_result_name = paste0("deconvoluted-", paths$cond1, ".txt"),
                        verbose = verbose,
                        pvalues_result_name = paste0("pvalues-", paths$cond1, ".txt"),
                        debug_seed = debug_seed,
                        threads = threads,
                        subsampling = subsampling,
                        subsampling_log = subsampling_log,
                        subsampling_num_pc = subsampling_num_pc,
                        subsampling_num_cells = subsampling_num_cells
    )
    run_cpdb_dot_plot(means_path = paste0(input_dir, "/cpdb_results_", paths$cond1, "/means-", paths$cond1, ".txt"),
                      pvalues_path = paste0(input_dir, "/cpdb_results_", paths$cond1, "/pvalues-", paths$cond1, ".txt"),
                      output_path = paste0(input_dir, "/cpdb_results_", paths$cond1),
                      output_name = paste0("dot_plot_", paths$cond1, ".pdf"),
                      rows = NULL,
                      columns = NULL,
                      verbose = verbose)
    run_cpdb_heatmap(metadata_path = paths$metadata_path1,
                     pvalues_path = paste0(input_dir, "/cpdb_results_", paths$cond1, "/pvalues-", paths$cond1, ".txt"),
                     output_path = paste0(input_dir, "/cpdb_results_", paths$cond1),
                     count_name = paste0("heatmap_count_", paths$cond1, ".pdf"),
                     log_name = paste0("heatmap_log_count_", paths$cond1, ".pdf"),
                     count_network_name = paste0("network_", paths$cond1, ".txt"),
                     interaction_count_name = paste0("interaction_count_", paths$cond1, ".txt"),
                     verbose = verbose)
    run_cpdb_from_files(data_path = paths$data_path2,
                        metadata_path = paths$metadata_path2,
                        method = method,
                        project_name = paste0("cpdb_results_", paths$cond2),
                        iterations = iterations,
                        threshold = threshold,
                        result_precision = result_precision,
                        counts_data = counts_data,
                        output_path = input_dir,
                        output_format = output_format,
                        means_result_name = paste0("means-", paths$cond2, ".txt"),
                        significant_means_result_name = paste0("significant-means-", paths$cond2, ".txt"),
                        deconvoluted_result_name = paste0("deconvoluted-", paths$cond2, ".txt"),
                        verbose = verbose,
                        pvalues_result_name = paste0("pvalues-", paths$cond2, ".txt"),
                        debug_seed = debug_seed,
                        threads = threads,
                        subsampling = subsampling,
                        subsampling_log = subsampling_log,
                        subsampling_num_pc = subsampling_num_pc,
                        subsampling_num_cells = subsampling_num_cells
    )
    run_cpdb_dot_plot(means_path = paste0(input_dir, "/cpdb_results_", paths$cond2, "/means-", paths$cond2, ".txt"),
                      pvalues_path = paste0(input_dir, "/cpdb_results_", paths$cond2, "/pvalues-", paths$cond2, ".txt"),
                      output_path = paste0(input_dir, "/cpdb_results_", paths$cond2),
                      output_name = paste0("dot_plot_", paths$cond2, ".pdf"),
                      rows = NULL,
                      columns = NULL,
                      verbose = verbose)
    run_cpdb_heatmap(metadata_path = paths$metadata_path2,
                     pvalues_path = paste0(input_dir, "/cpdb_results_", paths$cond2, "/pvalues-", paths$cond2, ".txt"),
                     output_path = paste0(input_dir, "/cpdb_results_", paths$cond2),
                     count_name = paste0("heatmap_count_", paths$cond2, ".pdf"),
                     log_name = paste0("heatmap_log_count_", paths$cond2, ".pdf"),
                     count_network_name = paste0("network_", paths$cond2, ".txt"),
                     interaction_count_name = paste0("interaction_count_", paths$cond2, ".txt"),
                     verbose = verbose)
  }
}


