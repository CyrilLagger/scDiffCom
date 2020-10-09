#' Run CellPhoneDB from a Seurat object
#'
#' @param seurat_object x
#' @param assay x
#' @param slot x
#' @param log_scale x
#' @param celltype_col_id x
#' @param input_species x
#' @param condition_col_id x
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
  celltype_col_id,
  input_species = "mouse",
  min_cells = 5,
  condition_col_id = NULL,
  input_dir = getwd(),
  create_plots = FALSE,
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
  subsampling_num_cells = NULL,
  return_full_dt = TRUE
) {
  message("Create file directory if not already existing.")
  if(!dir.exists(input_dir)) {
    dir.create(input_dir)
  }
  message("Create input files from Seurat to be used by CellPhoneDB.")
  paths <- create_cpdb_input(
    seurat_object = seurat_object,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    celltype_col_id = celltype_col_id,
    input_species = input_species,
    condition_col_id = condition_col_id,
    min_cells = min_cells,
    input_dir = input_dir
  )
  if(is.null(condition_col_id)) {
    n_run = 1
    project_name <- "cpdb_results_noCond"
    means_result_name <- NULL
    significant_means_result_name <- NULL
    econvoluted_result_name <-  NULL
    pvalues_result_name <- NULL
    message("Run CellphoneDB once (no condition).")
  } else {
    n_run = 2
    project_name <- paste0("cpdb_results_", paths$conds)
    means_result_name <- paste0("means-", paths$conds, ".txt")
    significant_means_result_name <- paste0("significant-means-", paths$conds, ".txt")
    deconvoluted_result_name <- paste0("deconvoluted-", paths$conds, ".txt")
    pvalues_result_name <- paste0("pvalues-", paths$conds, ".txt")
    message("Run CellPhoneDB on the two conditions.")
  }
  for(i in 1:n_run) {
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
  if(return_full_dt) {
    if(method != 'statistical_analysis') {
      message("Not possible to create the full data.table for the selected method.")
    } else {
      message("Create and write full data.table.")
      full_dt <- create_cpdb_cci(
        input_dir = input_dir,
        conds = paths$conds
        )
      if(is.null(condition_col_id)) {
        output_dir_full <-  paste0(input_dir, "/cpdb_full_table_noCond.txt")
      } else {
        output_dir_full <-  paste0(input_dir, "/cpdb_full_table_withCond.txt")
      }
      utils::write.table(
        full_dt,
        file = output_dir_full,
        quote = FALSE,
        col.names = TRUE,
        row.names = FALSE,
        sep = '\t'
      )
    }
  }
  if(input_species == "mouse") {
    message(paste0("Write human-mouse orthologs in ", input_dir, "/cpdb_human_mouse_orthologs.txt"))
    utils::write.table(
      paths$gene_mapping,
      file = paste0(input_dir, "/cpdb_human_mouse_orthologs.txt"),
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = '\t'
    )
  }
}

#' Extract a data matrix and a metadata data.frame from a Seurat object.
#'
#' @param seurat_object A Seurat object
#' @param assay The Seurat assay to pull data from; default is "RNA"
#' @param slot The Seurat slot to pull data from; default is "data"
#' @param log_scale Whether to return log-normalized or normalized data (only relevant when slot = "data"); default is TRUE
#' @param celltype_col_id Name of the column of the metadata data.frame containing the cell-type ids
#' @param input_species x
#' @param min_cells Minimum number of cells (per condition if relevant) required to keep a cell-type
#' @param condition_col_id Name of the column of the metadata data.frame containing the the condition on the cells. Set to NULL for no conditions
#' @param input_dir Directory path where to save the files that will be used as input for CellPhoneDB analysis
#'
#' @return Write the two files and return a list with the paths of the files and the names of the conditions (if relevant).
create_cpdb_input <- function(
  seurat_object,
  assay = "RNA",
  slot = "data",
  log_scale = FALSE,
  celltype_col_id,
  input_species = "mouse",
  min_cells = 5,
  condition_col_id = NULL,
  input_dir = getwd()
) {
  Gene <- cell_type <- NULL
  data <- extract_seurat_data(
    seurat_object = seurat_object,
    assay = assay,
    slot = slot,
    log_scale = log_scale,
    return_type = "data.table",
    verbose = TRUE
  )
  md <- extract_seurat_metadata(
    seurat_object = seurat_object,
    celltype_col_id = celltype_col_id,
    condition_col_id = condition_col_id
  )
  cell_type_filt <- filter_celltypes(
    metadata = md,
    min_cells = min_cells
    )
  md <- md[cell_type %in% cell_type_filt, ]
  cols <- colnames(data)[colnames(data) %in% c("rn", md$cell_id)]
  data <- data[, cols, with = FALSE]
  data.table::setnames(data, old = "rn", new = "Gene")
  data.table::setnames(md, old = c("cell_id"), new = c("Cell"))
  if(input_species == "mouse") {
    message("Converting mouse gene names to human gene names for CELLPHONEDB.")
    gene_mapping <- get_orthologs(
      genes = unique(data$Gene),
      input_species = "mouse",
      one2one = FALSE
    )
    data <- data.table::merge.data.table(
      data,
      gene_mapping[, c("mouse_symbol", "human_symbol")],
      by.x = "Gene",
      by.y = "mouse_symbol",
      all.x = TRUE
    )
    data <- stats::na.omit(data)
    data[, Gene := NULL]
    data.table::setnames(data, old = "human_symbol", new = "Gene")
    data.table::setcolorder(data, "Gene")
  } else if(input_species == "human") {
    message("Assuming genes symbols are human genes names.")
    gene_mapping <- NULL
  } else {
    stop("Species not supported: 'input_species' can be either 'human' or 'mouse'.")
  }
  if(is.null(condition_col_id)) {
    output_dir_data <- paste0(input_dir, "/cpdb_data_noCond.txt")
    output_dir_md <- paste0(input_dir, "/cpdb_metadata_noCond.txt")
    message(paste0("Writing CellPhoneDB input data to ", output_dir_data))
    utils::write.table(
      data,
      file = output_dir_data,
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = '\t'
    )
    message(paste0("Writing CellPhoneDB input metadata to ", output_dir_md))
    utils::write.table(
      md,
      file = output_dir_md,
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = '\t'
    )
    return(list(data_path = output_dir_data,
                md_path = output_dir_md,
                gene_mapping = gene_mapping))
  } else {
    conds <- unique(md$condition)
    if(length(conds) != 2) stop("Wrong number of groups in cell-type conditions (expected 2).")
    md1 <- md[md$condition == conds[[1]], ]
    md1$condition <- NULL
    md2 <- md[md$condition == conds[[2]], ]
    md2$condition <- NULL
    cols_1 <- colnames(data)[colnames(data) %in% c("Gene", md1$Cell)]
    cols_2 <- colnames(data)[colnames(data) %in% c("Gene", md2$Cell)]
    data1 <- data[, cols_1, with = FALSE]
    data2 <- data[, cols_2, with = FALSE]
    output_dir_data1 <- paste0(input_dir, "/cpdb_data_" , conds[[1]], ".txt")
    output_dir_md1 <- paste0(input_dir, "/cpdb_metadata_", conds[[1]], ".txt")
    output_dir_data2 <- paste0(input_dir, "/cpdb_data_" , conds[[2]], ".txt")
    output_dir_md2 <- paste0(input_dir, "/cpdb_metadata_", conds[[2]], ".txt")
    message(paste0("Writing CellPhoneDB input data to ", output_dir_data1))
    utils::write.table(
      data1,
      file = output_dir_data1,
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = '\t'
    )
    message(paste0("Writing CellPhoneDB input data to ", output_dir_data2))
    utils::write.table(
      data2,
      file = output_dir_data2,
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = '\t'
    )
    message(paste0("Writing CellPhoneDB input metadata to ", output_dir_md1))
    utils::write.table(
      md1,
      file = output_dir_md1,
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = '\t'
    )
    message(paste0("Writing CellPhoneDB input metadata to ", output_dir_md2))
    utils::write.table(
      md2,
      file = output_dir_md2,
      quote = FALSE,
      col.names = TRUE,
      row.names = FALSE,
      sep = '\t'
    )
    return(list(conds = c(conds[[1]], conds[[2]]),
                data_path = c(output_dir_data1, output_dir_data2),
                md_path = c(output_dir_md1, output_dir_md2),
                gene_mapping = gene_mapping))
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
run_cpdb_from_files <- function(
  data_path,
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
#' @param conds x
#'
#' @return x
create_cpdb_cci <- function(
  input_dir,
  conds = NULL
) {
  if(is.null(conds)) {
    path_res <- list(
      paste0(input_dir, '/cpdb_results_noCond/means.txt'),
      paste0(input_dir, '/cpdb_results_noCond/pvalues.txt')#,
      #paste0(input_dir, "/cpdb_results_noCond/deconvoluted.txt")
    )
    names(path_res) <- c("means_noCond", "pvalues_noCond")#, "deconv_noCond")
  } else {
    path_res <- list(
      paste0(input_dir, "/cpdb_results_", conds[[1]], "/means-", conds[[1]], ".txt"),
      paste0(input_dir, "/cpdb_results_", conds[[2]], "/means-", conds[[2]], ".txt"),
      paste0(input_dir, "/cpdb_results_", conds[[1]], "/pvalues-", conds[[1]], ".txt"),
      paste0(input_dir, "/cpdb_results_", conds[[2]], "/pvalues-", conds[[2]], ".txt")#,
     # paste0(input_dir, "/cpdb_results_", conds[[1]], "/deconvoluted-", conds[[1]], ".txt"),
      #paste0(input_dir, "/cpdb_results_", conds[[2]], "/deconvoluted-", conds[[2]], ".txt")
    )
    names(path_res) <- c(paste0("means_", conds), paste0("pvalues_", conds))#, paste0("deconv_", conds))
  }
  cols_to_rm <- c("interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b",
                  "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin")
  cpdb_res <- lapply(
    seq_along(path_res),
    function(i) {
      temp <- utils::read.table(
        path_res[[i]],
        header = TRUE,
        sep = "\t")
      data.table::setDT(temp)
      temp <- temp[, (cols_to_rm) := NULL]
      temp <- data.table::melt.data.table(
        temp,
        id.vars = "id_cp_interaction",
        variable.name = "cell_type_pair",
        value.name = names(path_res)[[i]]
      )
      return(temp)
    }
  )
  names(cpdb_res) <- names(path_res)
  if(is.null(conds)) {
    full_res <- data.table::merge.data.table(
      cpdb_res[[1]],
      cpdb_res[[2]],
      by = c("id_cp_interaction", "cell_type_pair")
    )
  } else {
    res_1 <- data.table::merge.data.table(
      cpdb_res[[1]],
      cpdb_res[[3]],
      by = c("id_cp_interaction", "cell_type_pair")
    )
    res_2 <- data.table::merge.data.table(
      cpdb_res[[2]],
      cpdb_res[[4]],
      by = c("id_cp_interaction", "cell_type_pair")
    )
    full_res <- data.table::merge.data.table(
      res_1,
      res_2,
      by = c("id_cp_interaction", "cell_type_pair"),
      all = TRUE
    )
  }
  return(full_res)
}

