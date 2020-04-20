#' Select LR genes present in the dataset
#'
#' @param data A matrix or data.frame with genes as rows and cells as columns
#' @param LR_data A data.frame with one LR pair per row
#'
#' @return A list containing the subsetted data and the subsetted LR-pairs
#' @export
#'
#' @examples
preprocess_LR <- function(
  data,
  LR_data
) {
  message(paste0("Number of considered LR pairs: ", length(unique(LR_data$SYMB_LR)), "."))
  LR_keep <- LR_data[LR_data$GENESYMB_L %in% rownames(data) &
                       LR_data$GENESYMB_R %in% rownames(data), ]
  LR_genes <- unique(c(unique(LR_keep$GENESYMB_L), unique(LR_keep$GENESYMB_R)))
  data_keep <- data[rownames(data) %in% LR_genes, ]
  message(paste0("Number of LR pairs in the dataset: ", length(unique(LR_keep$SYMB_LR)), "."))
  #stop("Transform LR data in a universal dataframe")
  LR_keep$ligand <- LR_keep$GENESYMB_L
  LR_keep$receptor <- LR_keep$GENESYMB_R
  return(list(data = data_keep, LR_df = LR_keep))
}
