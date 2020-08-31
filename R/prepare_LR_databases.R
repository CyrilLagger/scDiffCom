#' Create a dataframe with LR pairs from our 4 sources
#'
#' @param one2one logical indicating the orthology relationship used during conversion
#'
#' @return data.frame with LR pairs and to which sources they belong
#' @export
#'
#' @examples
#' \dontrun{
#' LRall_one2one <- aggregate_LR(one2one = TRUE)
#' LRall_many2many <- aggregate_LR(one2one = FALSE)
#' }
aggregate_LR <- function(
  one2one = FALSE
) {
  CONF_L <- CONF_L_nichenet <- CONF_L_scsr <-
    CONF_R <- CONF_R_nichenet <- CONF_R_scsr <-
    CONF_L.x <- CONF_R.x <- CONF_L.y <- CONF_R.y <-
    TYPE_L <- TYPE_L_nichenet <- TYPE_L_scsr <-
    TYPE_R <- TYPE_R_nichenet <- TYPE_R_scsr <-
    TYPE_L.x <- TYPE_R.x <- TYPE_L.y <- TYPE_R.y <-
    sctensor <- cpdb <- nichenet <- scsr <-
    source_cpdb <- SYMB_ab <- SYMB_ba <-  NULL
  #load each data.table
  scT_dt <- prepare_LR_scTensor(
    detailed = FALSE
  )
  scT_dt$sctensor <- TRUE
  scsr_dt <- prepare_LR_scsr(
    detailed = FALSE,
    one2one = one2one
  )
  scsr_dt$scsr <- TRUE
  niche_dt <- prepare_LR_nichenet(
    detailed = FALSE,
    one2one = one2one
  )
  niche_dt$nichenet <- TRUE
  cpdb_dt <- prepare_LR_cpdb(
    one2one = one2one,
    deconvoluted = TRUE
  )
  cpdb_dt$cpdb <- TRUE
  #merge the LR pairs
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
    sort = FALSE,
    suffixes = c("_scsr", "_nichenet")
  )
  dt <- data.table::merge.data.table(
    dt,
    cpdb_dt[,c("SYMB_ab", "cpdb")],
    by.x = c("SYMB_LR"),
    by.y = c( "SYMB_ab"),
    all.x = TRUE,
    sort = FALSE
  )
  dt <- data.table::merge.data.table(
    dt,
    cpdb_dt[,c("SYMB_ba", "cpdb")],
    by.x = c("SYMB_LR"),
    by.y = c("SYMB_ba"),
    all.x = TRUE,
    sort = FALSE
  )
  dt$cpdb <- dt$cpdb.x | dt$cpdb.y
  dt[,c("cpdb.x", "cpdb.y") := NULL]
  dt[, CONF_L := ifelse(
    !is.na(CONF_L_nichenet),
    CONF_L_nichenet,
    ifelse(
      !is.na(CONF_L_scsr),
      CONF_L_scsr,
      ifelse(
        sctensor == TRUE,
        1,
        NA
      )
    )) ]
  dt[, CONF_R := ifelse(
    !is.na(CONF_R_nichenet),
    CONF_R_nichenet,
    ifelse(
      !is.na(CONF_R_scsr),
      CONF_R_scsr,
      ifelse(
        sctensor == TRUE,
        1,
        NA
      )
    )) ]
  dt[, TYPE_L := ifelse(
    !is.na(TYPE_L_nichenet),
    TYPE_L_nichenet,
    ifelse(
      !is.na(TYPE_L_scsr),
      TYPE_L_scsr,
      ifelse(
        sctensor == TRUE,
        "sctensor",
        NA
      )
    )) ]
  dt[, TYPE_R := ifelse(
    !is.na(TYPE_R_nichenet),
    TYPE_R_nichenet,
    ifelse(
      !is.na(TYPE_R_scsr),
      TYPE_R_scsr,
      ifelse(
        sctensor == TRUE,
        "sctensor",
        NA
      )
    )) ]
  dt[,c("CONF_L_nichenet", "CONF_R_nichenet", "CONF_L_scsr", "CONF_R_scsr",
        "TYPE_L_nichenet", "TYPE_R_nichenet", "TYPE_L_scsr", "TYPE_R_scsr") := NULL]
  #There is some uncertainty in the ordering of cpdb LR pairs that are not present in the rest of the data
  #We order the ones that are inconsistent with the rest of the data
  all_genes <- unique(c(dt$GENESYMB_L, dt$GENESYMB_R))
  onlyL_genes <- setdiff(
    unique(dt$GENESYMB_L),
    unique(dt$GENESYMB_R)
  )
  onlyR_genes <- setdiff(
    unique(dt$GENESYMB_R),
    unique(dt$GENESYMB_L)
  )
  common_genes <- intersect(
    unique(dt$GENESYMB_L),
    unique(dt$GENESYMB_R)
  )
  cpdb_only_dt <- cpdb_dt[!(SYMB_ab %in% dt$SYMB_LR) & !(SYMB_ba %in% dt$SYMB_LR), ]
  cond <- cpdb_only_dt$GENESYMB_a %in% onlyR_genes & !(cpdb_only_dt$GENESYMB_b %in% onlyR_genes) |
               cpdb_only_dt$GENESYMB_b %in% onlyL_genes & !(cpdb_only_dt$GENESYMB_a %in% onlyL_genes)
  cpdb_only_dt$SYMB_LR <- ifelse(
    cond,
    cpdb_only_dt$SYMB_ba,
    cpdb_only_dt$SYMB_ab
  )
  cpdb_only_dt$GENESYMB_L <- ifelse(
    cond,
    cpdb_only_dt$GENESYMB_b,
    cpdb_only_dt$GENESYMB_a
  )
  cpdb_only_dt$GENESYMB_R <- ifelse(
    cond,
    cpdb_only_dt$GENESYMB_a,
    cpdb_only_dt$GENESYMB_b
  )
  cpdb_only_dt$CONF_L <- ifelse(
    cond,
    cpdb_only_dt$CONF_b,
    cpdb_only_dt$CONF_a
  )
  cpdb_only_dt$CONF_R <- ifelse(
    cond,
    cpdb_only_dt$CONF_a,
    cpdb_only_dt$CONF_b
  )
  cpdb_only_dt$TYPE_L <- ifelse(
    cond,
    cpdb_only_dt$TYPE_b,
    cpdb_only_dt$TYPE_a
  )
  cpdb_only_dt$TYPE_R <- ifelse(
    cond,
    cpdb_only_dt$TYPE_a,
    cpdb_only_dt$TYPE_b
  )
  dt <- data.table::merge.data.table(
    dt,
    cpdb_only_dt[,c("GENESYMB_L", "GENESYMB_R","SYMB_LR","cpdb", "CONF_L", "CONF_R", "TYPE_L", "TYPE_R")],
    by= c("GENESYMB_L", "GENESYMB_R", "SYMB_LR"),
    all = TRUE,
    sort = FALSE
  )
  dt$cpdb <- dt$cpdb.x | dt$cpdb.y
  dt[, c("cpdb.x", "cpdb.y") := NULL]
  dt[, CONF_L := ifelse(
    !is.na(CONF_L.x),
    CONF_L.x,
    CONF_L.y
  )]
  dt[, CONF_R := ifelse(
    !is.na(CONF_R.x),
    CONF_R.x,
    CONF_R.y
  )]
  dt[, TYPE_L := ifelse(
    !is.na(TYPE_L.x),
    TYPE_L.x,
    TYPE_L.y
  )]
  dt[, TYPE_R := ifelse(
    !is.na(TYPE_R.x),
    TYPE_R.x,
    TYPE_R.y
  )]
  dt[,c("CONF_L.x", "CONF_R.x", "CONF_L.y", "CONF_R.y",
        "TYPE_L.x", "TYPE_R.x", "TYPE_L.y", "TYPE_R.y") := NULL]
  dt[, source_cpdb := ifelse(cpdb == TRUE, "cpdb", NA)]
  dt[, scsr := ifelse(is.na(scsr), FALSE, TRUE)]
  dt[, nichenet := ifelse(is.na(nichenet), FALSE, TRUE)]
  dt[, cpdb := ifelse(is.na(cpdb), FALSE, TRUE)]
  dt[, sctensor := ifelse(is.na(sctensor), FALSE, TRUE)]
  data.table::setcolorder(
    x = dt,
    neworder = c("GENESYMB_L", "GENESYMB_R", "SYMB_LR",
                 "scsr", "cpdb", "nichenet", "sctensor",
                 "source_scsr", "source_cpdb", "source_nichenet", "source_sctensor",
                 "CONF_L", "TYPE_L", "CONF_R", "TYPE_R")
  )
  return(dt)
}

#' Create a LR data.table based on the scTensor mouse database
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
  LR_L_match <- AnnotationDbi::select(
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
    LR$LIGAND <- LR_L_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  if(identical(LR_R_match$ENTREZID, as.character(LR$GENEID_R))) {
    LR$RECEPTOR <- LR_R_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  data.table::setDT(LR)
  LR$LR_ID <- paste(LR$LIGAND, LR$RECEPTOR, sep = "_")
  if(!detailed) {
    LR <- LR[, c("LR_ID", "LIGAND", "RECEPTOR", "SOURCEDB")]
    LR <- LR[!duplicated(LR$LR_ID), ]
    data.table::setnames(LR, old = "SOURCEDB", new = "SOURCE_SCTENSOR")
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
  one2one = FALSE
) {
  LR <- SingleCellSignalR::LRdb
  data.table::setDT(LR)
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
  LR <- data.table::merge.data.table(
    LR,
    L,
    by.x = "ligand",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  LR <- data.table::merge.data.table(
    LR,
    R,
    by.x = "receptor",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setnames(
    x = LR,
    old = c("mouse_symbol.x", "confidence.x", "type.x", "mouse_symbol.y", "confidence.y", "type.y"),
    new = c("LIGAND", "CONF_L", "TYPE_L", "RECEPTOR", "CONF_R", "TYPE_R")
  )
  if(!detailed) {
    LR <- stats::na.omit(LR)
    LR <- LR[, c("LIGAND", "RECEPTOR", "source", "CONF_L", "TYPE_L", "CONF_R", "TYPE_R")]
    LR$LR_ID <- paste(LR$LIGAND, LR$RECEPTOR, sep = "_")
    LR <- LR[!duplicated(LR$LR_ID), ]
    data.table::setnames(x = LR, old = "source", new =  "source_scsr")
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
  data.table::setDT(niche)
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
  niche <- data.table::merge.data.table(
    niche,
    L,
    by.x = "from",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  niche <- data.table::merge.data.table(
    niche,
    R,
    by.x = "to",
    by.y = "human_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setnames(
    x = niche,
    old = c("mouse_symbol.x", "confidence.x", "type.x", "mouse_symbol.y", "confidence.y", "type.y"),
    new = c("GENESYMB_L", "CONF_L", "TYPE_L", "GENESYMB_R", "CONF_R", "TYPE_R")
  )
  if(!detailed) {
    niche <- stats::na.omit(niche)
    niche <- niche[, c("GENESYMB_L", "GENESYMB_R", "database", "CONF_L", "TYPE_L", "CONF_R", "TYPE_R")]
    niche$SYMB_LR <- paste(niche$GENESYMB_L, niche$GENESYMB_R, sep = "_")
    niche <- niche[!duplicated(niche$SYMB_LR), ]
    setnames(niche, "database", "source_nichenet")
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
  data.table::setDT(LR)
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
    LR <- data.table::merge.data.table(
      LR,
      Ga,
      by.x = "SYMB_a",
      by.y = "human_symbol",
      all.x = TRUE,
      sort = FALSE
    )
    LR <- data.table::merge.data.table(
      LR,
      Gb,
      by.x = "SYMB_b",
      by.y = "human_symbol",
      all.x = TRUE,
      sort = FALSE
    )
    data.table::setnames(
      x = LR,
      old = c("mouse_symbol.x", "confidence.x", "type.x", "mouse_symbol.y", "confidence.y", "type.y"),
      new = c("GENESYMB_a", "CONF_a", "TYPE_a", "GENESYMB_b", "CONF_b", "TYPE_b")
    )
    LR <- stats::na.omit(LR)
    LR <- LR[, c("GENESYMB_a", "GENESYMB_b", "CONF_a", "TYPE_a", "CONF_b", "TYPE_b")]
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

