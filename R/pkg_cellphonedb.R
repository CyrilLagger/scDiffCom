#' Title
#'
#' @param seurat_obj a Seurat object
#' @param assay assay to pull data from
#' @param slot slot to pull data from
#' @param log_scale logical
#' @param seurat_cell_type_id xx
#' @param condition_id xx
#' @param output_dir xx
#'
#' @return
#' @export
#'
#' @examples
create_cpdb_input <- function(seurat_obj,
                                  assay = "RNA",
                                  slot = "data",
                                  log_scale = TRUE,
                                  seurat_cell_type_id,
                                  condition_id = NULL,
                                  output_dir
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
    message(paste0("Writing CellPhoneDB input data to ", output_dir, "/cpdb_data_noCond.txt"))
    utils::write.table(data,
                       file = paste0(output_dir, "/cpdb_data_noCond.txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')
    message(paste0("Writing CellPhoneDB input metadata to ", output_dir, "/cpdb_metadata_noCond.txt"))
    utils::write.table(metadata,
                       file = paste0(output_dir, "/cpdb_metadata_noCond.txt"),
                       quote = FALSE,
                       col.names = TRUE,
                       row.names = FALSE,
                       sep = '\t')
  } else {
    stop("Applying CellPhoneDB to multiple conditions is not yet implemented. Stay tuned!")
  }
}


# run_cpdb_from_files <- function(data_path,
#                                 metadata_path
# ) {
#   #system('cellphonedb method statistical_analysis metadata.txt counts.txt --iterations=10 --threads=2')
#
#   #system('cellphonedb plot dot_plot')
#
#   #system('cellphonedb plot heatmap_plot metadata.txt')
# }
