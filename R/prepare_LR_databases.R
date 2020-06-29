#' Create a dataframe with LR pairs from our 4 sources
#'
#' @param one2one logical indicating the orthology relationship used during conversion
#' @param cpdb_usage logical indicating if returning cpdb
#'
#' @return data.frame with LR pairs and to which sources they belong
#' @export
#'
#' @examples
#' \dontrun{
#' LRall_one2one <- aggregate_LR(one2one = TRUE, cpdb = TRUE)
#' LRall_many2many <- aggregate_LR(one2one = FALSE, cpdb = TRUE)
#' }
aggregate_LR <- function(
  one2one = TRUE,
  cpdb_usage = TRUE
) {
  scT_dt <- prepare_LR_scTensor(
    detailed = FALSE
  )
  scT_dt$sctensor <- TRUE
  data.table::setDT(scT_dt)
  scsr_dt <- prepare_LR_scsr(
    detailed = FALSE,
    one2one = one2one
  )
  scsr_dt$scsr <- TRUE
  data.table::setDT(scsr_dt)
  niche_dt <- prepare_LR_nichenet(
    detailed = FALSE,
    one2one = one2one
  )
  niche_dt$nichenet <- TRUE
  data.table::setDT(niche_dt)
  dt <- data.table::merge.data.table(
    scT_dt,
    scsr_dt,
    by = c("GENESYMB_L", "GENESYMB_R", "SYMB_LR"),
    all = TRUE,
    sort = FALSE
  )
  dt <- data.table::merge.data.table(
    dt,
    niche_dt,
    by = c("GENESYMB_L", "GENESYMB_R", "SYMB_LR"),
    all = TRUE,
    sort = FALSE
  )
  dt[is.na(dt)] <- FALSE
  if(cpdb_usage) {
    cpdb_dt <- prepare_LR_cpdb(
      one2one = one2one,
      deconvoluted = TRUE
    )
    cpdb_dt$cpdb <- TRUE
    data.table::setDT(cpdb_dt)
    dt <- data.table::merge.data.table(
      dt,
      cpdb_dt[,c("SYMB_ab", "source_cpdb", "cpdb")],
      by.x = c("SYMB_LR"),
      by.y = c( "SYMB_ab"),
      all.x = TRUE,
      sort = FALSE
    )
    dt[is.na(dt)] <- FALSE
    dt <- data.table::merge.data.table(
      dt,
      cpdb_dt[,c("SYMB_ba", "cpdb")],
      by.x = c("SYMB_LR"),
      by.y = c("SYMB_ba"),
      all.x = TRUE,
      sort = FALSE
    )
    dt[is.na(dt)] <- FALSE
    dt$cpdb <- dt$cpdb.x | dt$cpdb.y
    dt[,c("cpdb.x", "cpdb.y") := NULL]
    not_cpdb_dt <- cpdb_dt[!(cpdb_dt$SYMB_ab %in% dt$SYMB_LR) & !(cpdb_dt$SYMB_ba %in% dt$SYMB_LR), ]
    dt <- data.table::merge.data.table(
      dt,
      not_cpdb_dt[,c("SYMB_ab","cpdb")],
      by.x = c("SYMB_LR"),
      by.y = c("SYMB_ab"),
      all = TRUE,
      sort = FALSE
    )
    dt[is.na(dt)] <- FALSE
    dt$cpdb <- dt$cpdb.x | dt$cpdb.y
    dt[,c("cpdb.x", "cpdb.y") := NULL]
  }
  dt$source_cpdb <- ifelse(dt$cpdb == TRUE, "cpdb", FALSE )
  return(dt)
}

#' Create a LR data.frame based on the scTensor database
#'
#' @param detailed logical indicating if returning only 3 columns or the more detailed data.frame.
#'
#' @return data.frame with 3 or 6 columns
#' @export
#'
#' @examples
#' \dontrun{
#' LRsct <- prepare_LR_scTensor(detailed = FALSE)
#' }
prepare_LR_scTensor <- function(
  detailed = FALSE
) {
  key <- AnnotationDbi::keys(
    LRBase.Mmu.eg.db::LRBase.Mmu.eg.db,
    keytype="GENEID_L"
  )
  LR <- AnnotationDbi::select(
    LRBase.Mmu.eg.db::LRBase.Mmu.eg.db,
    keys = key,
    columns = c("GENEID_L", "GENEID_R", "SOURCEDB"),
    keytype = "GENEID_L"
  )
  LR <- unique(LR)
  #change the names of the genes from EntrezID to Symbol
  ah <- AnnotationHub::AnnotationHub()
  hs <- AnnotationHub::query(
    ah,
    c("OrgDb", "Mus musculus")
  )[[1]]

  LR_L_match <-AnnotationDbi::select(
    hs,
    column=c("SYMBOL", "ENTREZID"),
    keytype="ENTREZID",
    keys= as.character(LR$GENEID_L)
  )
  LR_R_match <- AnnotationDbi::select(
    hs,
    column=c("SYMBOL", "ENTREZID"),
    keytype="ENTREZID",
    keys= as.character(LR$GENEID_R)
  )
  if(identical(LR_L_match$ENTREZID, as.character(LR$GENEID_L))) {
    LR$GENESYMB_L <- LR_L_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  if(identical(LR_R_match$ENTREZID, as.character(LR$GENEID_R))) {
    LR$GENESYMB_R <- LR_R_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  LR$SYMB_LR <- paste(LR$GENESYMB_L, LR$GENESYMB_R, sep = "_")
  if(!detailed) {
    LR <- LR[, c("GENESYMB_L", "GENESYMB_R", "SYMB_LR", "SOURCEDB")]
    LR <- LR[!duplicated(LR$SYMB_LR), ]
    colnames(LR)[colnames(LR) == "SOURCEDB"] <- "source_sctensor"
  }
  return(LR)
}

#' Create a LR data.frame based on SingleCellSignalR database
#'
#' @param detailed logical indicating if returning only 3 columns or the more detailed data.frame.
#' @param one2one logical indicating if using one2one orthology relationship
#'
#' @return data.frame with 3 LR columns and possibly more detailed information
#' @export
#'
#' @examples
#' \dontrun{
#' LRdb_one2one <- prepare_LR_scsr(detailed = FALSE, one2one = TRUE)
#' LRdb_many2many <- prepare_LR_scsr(detailed = FALSE, one2one = FALSE)
#' }
prepare_LR_scsr <- function(
  detailed = FALSE,
  one2one = TRUE
) {
  LR <- SingleCellSignalR::LRdb
  L <- get_orthologs(
    unique(LR$ligand),
    input_species = "human",
    one2one = one2one
  )
  L <- stats::na.omit(L)
  R <- get_orthologs(
    unique(LR$receptor),
    input_species = "human",
    one2one = one2one
  )
  R <- stats::na.omit(R)
  LR <- merge(
    LR,
    L,
    by.x = "ligand",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  LR <- merge(
    LR,
    R,
    by.x = "receptor",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  colnames(LR)[14:15] <- c("GENESYMB_L", "GENESYMB_R")
  if(!detailed) {
    LR <- stats::na.omit(LR)
    LR <- LR[, c("GENESYMB_L", "GENESYMB_R", "source")]
    LR$SYMB_LR <- paste(LR$GENESYMB_L, LR$GENESYMB_R, sep = "_")
    LR <- LR[!duplicated(LR$SYMB_LR), ]
    colnames(LR)[colnames(LR) == "source"] <- "source_scsr"
  }
  return(LR)
}

#' Create a LR data.frame based on nicheNet database
#'
#' @param detailed logical indicating if returning only 3 columns or the more detailed data.frame.
#' @param one2one logical indicating if using one2one orthology relationship
#'
#' @return data.frame with 3 LR columns and possibly more detailed information
#' @export
#'
#' @examples
#' \dontrun{
#' LRniche_one2one <- prepare_LR_nichenet(detailed = FALSE, one2one = TRUE)
#' LRniche_many2many <- prepare_LR_nichenet(detailed = FALSE, one2one = FALSE)
#' }
prepare_LR_nichenet <- function(
  detailed = FALSE,
  one2one = TRUE
) {
  niche <- nichenetr::lr_network
  L <- get_orthologs(
    genes = unique(niche$from),
    input_species = "human",
    one2one = one2one
  )
  L <- stats::na.omit(L)
  R <- get_orthologs(
    genes = unique(niche$to),
    input_species = "human",
    one2one = one2one
  )
  R <- stats::na.omit(R)
  niche <- merge(
    niche,
    L,
    by.x = "from",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  niche <- merge(
    niche,
    R,
    by.x = "to",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  colnames(niche)[5:6] <- c("GENESYMB_L", "GENESYMB_R")
  if(!detailed) {
    niche <- stats::na.omit(niche)
    niche <- niche[, c("GENESYMB_L", "GENESYMB_R", "database")]
    niche$SYMB_LR <- paste(niche$GENESYMB_L, niche$GENESYMB_R, sep = "_")
    niche <- niche[!duplicated(niche$SYMB_LR), ]
    colnames(niche)[colnames(niche) == "database"] <- "source_nichenet"
  }
  return(niche)
}

#' Create a LR data.frame based on CellPhoneDB database.
#'
#' @param one2one logical indicating if using one2one orthology relationship
#' @param deconvoluted logical indicating if using 1:1 LR versus complex LR
#'
#' @return data.frame with 3 LR columns and possibly more detailed information
#' @export
#'
#' @examples
#' \dontrun{
#' LRcpdb_one2one <- prepare_LR_cpdb(one2one = TRUE, deconvoluted = TRUE)
#' LRcpdb_many2many <- prepare_LR_cpdb(one2one = FALSE, deconvoluted = TRUE)
#' }
prepare_LR_cpdb <- function(
  one2one = TRUE,
  deconvoluted = TRUE
) {
  LR <- create_LR_cpdb(
    deconvoluted = deconvoluted
  )
  if(deconvoluted) {
    Ga <- get_orthologs(
      genes = unique(LR$SYMB_a),
      input_species = "human",
      one2one = one2one
    )
    Ga <- stats::na.omit(Ga)
    Gb <- get_orthologs(
      genes = unique(LR$SYMB_b),
      input_species = "human",
      one2one = one2one
    )
    Gb <- stats::na.omit(Gb)
    LR <- merge(
      LR,
      Ga,
      by.x = "SYMB_a",
      by.y = "human_symbol",
      all.x = TRUE,
      sort = FALSE
    )
    LR <- merge(
      LR,
      Gb,
      by.x = "SYMB_b",
      by.y = "human_symbol",
      all.x = TRUE,
      sort = FALSE
    )
    colnames(LR)[4:5] <- c("GENESYMB_a", "GENESYMB_b")
    LR <- stats::na.omit(LR)
    LR <- LR[, c(4:5)]
    LR$SYMB_ab <- paste(LR$GENESYMB_a, LR$GENESYMB_b, sep = "_")
    LR <- LR[!duplicated(LR$SYMB_ab), ]
    LR$SYMB_ba <- paste(LR$GENESYMB_b, LR$GENESYMB_a, sep = "_")
  } else {
    stop("Not supported yet.")
    La <- get_orthologs(
      genes = stats::na.omit(unique(LR$gene_1_a)),
      input_species = "human",
      one2one = one2one
    )
    La <- stats::na.omit(La)
    Lb <- get_orthologs(
      genes = stats::na.omit(unique(LR$gene_1_b)),
      input_species = "human",
      one2one = one2one
    )
    Lb <- stats::na.omit(Lb)
    Lc <- get_orthologs(
      genes = stats::na.omit(unique(LR$gene_1_c)),
      input_species = "human",
      one2one = one2one
    )
    Lc <- stats::na.omit(Lc)
    Ra <- get_orthologs(
      genes = stats::na.omit(unique(LR$gene_2_a)),
      input_species = "human",
      one2one = one2one
    )
    Ra <- stats::na.omit(Ra)
    Rb <- get_orthologs(
      genes = stats::na.omit(unique(LR$gene_2_b)),
      input_species = "human",
      one2one = one2one
    )
    Rb <- stats::na.omit(Rb)
    Rc <- get_orthologs(
      genes = stats::na.omit(unique(LR$gene_2_c)),
      input_species = "human",
      one2one = one2one
    )
    Rc <- stats::na.omit(Rc)
  }
  LR$source_cpdb <- "cpdb"
  return(LR)
}

#' Convert CellphoneDB database in a dataframe either complex or deconvoluted
#'
#' @param deconvoluted logical indicating if returning the complex LR with 1 to 3 genes per category or the 1:1 deconvoluted data.frame.
#'
#' @return data.frame
create_LR_cpdb <- function(
  deconvoluted = TRUE
) {
  data <- LRcp_raw
  res <- sapply(unique(data$interaction_table$id_cp_interaction), function(id_cp) {
    multi_id_1 <- data$interaction_table[data$interaction_table$id_cp_interaction == id_cp, "multidata_1_id"]
    multi_id_2 <- data$interaction_table[data$interaction_table$id_cp_interaction == id_cp, "multidata_2_id"]

    if(data$multidata_table[data$multidata_table$id_multidata == multi_id_1, "is_complex"] == 0) {
      id_prot_1 <- data$protein_table[data$protein_table$protein_multidata_id == multi_id_1, "id_protein"]
      gene_1_a <- data$gene_table[data$gene_table$protein_id == id_prot_1, "hgnc_symbol"][[1]]
      gene_1_b <- NA
      gene_1_c <- NA
    } else {
      comp_df_1 <- data$complex_composition_table[data$complex_composition_table$complex_multidata_id == multi_id_1,]
      gene_1_a <- data$gene_table[data$gene_table$protein_id == comp_df_1[1,"protein_multidata_id"], "hgnc_symbol"][[1]]
      gene_1_b <- data$gene_table[data$gene_table$protein_id == comp_df_1[2,"protein_multidata_id"], "hgnc_symbol"][[1]]
      if(nrow(comp_df_1) == 3) {
        gene_1_c <- data$gene_table[data$gene_table$protein_id == comp_df_1[3,"protein_multidata_id"], "hgnc_symbol"]
      } else {
        gene_1_c <- NA
      }

    }
    if(data$multidata_table[data$multidata_table$id_multidata == multi_id_2, "is_complex"] == 0) {
      id_prot_2 <- data$protein_table[data$protein_table$protein_multidata_id == multi_id_2, "id_protein"]
      gene_2_a <- data$gene_table[data$gene_table$protein_id == id_prot_2, "hgnc_symbol"][[1]]
      gene_2_b <- NA
      gene_2_c <- NA
    } else {
      comp_df_2 <- data$complex_composition_table[data$complex_composition_table$complex_multidata_id == multi_id_2,]
      gene_2_a <- data$gene_table[data$gene_table$protein_id == comp_df_2[1,"protein_multidata_id"], "hgnc_symbol"][[1]]
      gene_2_b <- data$gene_table[data$gene_table$protein_id == comp_df_2[2,"protein_multidata_id"], "hgnc_symbol"][[1]]
      if(nrow(comp_df_2) == 3) {
        gene_2_c <- data$gene_table[data$gene_table$protein_id == comp_df_2[3,"protein_multidata_id"], "hgnc_symbol"]
      } else {
        gene_2_c <- NA
      }

    }
    return(c("gene_1_a" = gene_1_a,
             "gene_1_b" = gene_1_b,
             "gene_1_c" = gene_1_c,
             "gene_2_a" = gene_2_a,
             "gene_2_b" = gene_2_b,
             "gene_2_c" = gene_2_c))

  })
  res <- as.data.frame(t(res), stringsAsFactors = FALSE)
  if(deconvoluted) {
    res <- do.call("rbind", lapply(1:3, function(i) {
      do.call("rbind",lapply(4:6, function(j) {
        df <- stats::na.omit(res[,c(i,j)])
        colnames(df) <- c("SYMB_a", "SYMB_b")
        df$LR <- paste(df$SYMB_a, df$SYMB_b, sep = "_")
        df <- df[!duplicated(df$LR),]
      }) )
    }))
    res <- res[!duplicated(res$LR),]
  }
  return(res)
}

