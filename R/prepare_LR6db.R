#' Create a data.table with the ligand-receptor interactions from 6 databases.
#'
#' @param one2one logical indicating if the orthology conversion has to be limited to one2one homology; default is FALSE.
#' @param curated logical indicating if only the curated LR interactions should be returned.
#'
#' @return data.table with ligand-receptor interactions and relevant properties from 6 databases.
combine_LR_db <- function(
  one2one = FALSE,
  curated = TRUE
) {
  DATABASE <- SOURCE <- ANNOTATION <- FAMILY <- SUBFAMILY <- keep_subLR <- NULL
  LR_sct <- prepare_LR_scTensor()
  LR_scsr <- prepare_LR_scsr(one2one = one2one)
  #LR_niche <- prepare_LR_nichenet(one2one = one2one)
  LR_cpdb <- prepare_LR_cpdb(one2one = one2one, deconvoluted = FALSE)
  #LR_cc <- prepare_LR_CellChat()
  LR_ic <- prepare_LR_ICELLNET(one2one = one2one)
  LR_full <- rbindlist(
    list(
      "SCTENSOR" = LR_sct,
      "SCSR" = LR_scsr,
      "NICHENET" = LR_niche,
      "CELLPHONEDB" = LR_cpdb,
      "CELLCHAT" = LR_cc,
      "ICELLNET" = LR_ic),
    use.names = TRUE,
    fill = TRUE,
    idcol = "DATABASE"
  )
  col_md <- colnames(LR_full)
  col_md <- col_md[col_md != "LR_SORTED"]
  db_names <- c("CELLCHAT", "CELLPHONEDB", "SCSR", "SCTENSOR", "ICELLNET", "NICHENET")
  LR_full <- data.table::dcast.data.table(
    LR_full,
    formula = LR_SORTED ~ DATABASE,
    value.var = col_md
  )
  LR_full[, paste0("LIGAND_", 1:2) := lapply(
    1:2,
    function(i) {
      ifelse(
        !is.na(get(paste0("LIGAND_", i, "_CELLPHONEDB"))),
        get(paste0("LIGAND_", i, "_CELLPHONEDB")),
        ifelse(
          !is.na(get(paste0("LIGAND_", i, "_ICELLNET"))),
          get(paste0("LIGAND_", i, "_ICELLNET")),
          ifelse(
            !is.na(get(paste0("LIGAND_", i, "_SCSR"))),
            get(paste0("LIGAND_", i, "_SCSR")),
            ifelse(
              !is.na(get(paste0("LIGAND_", i, "_SCTENSOR"))),
              get(paste0("LIGAND_", i, "_SCTENSOR")),
              ifelse(
                !is.na(get(paste0("LIGAND_", i, "_NICHENET"))),
                get(paste0("LIGAND_", i, "_NICHENET")),
                get(paste0("LIGAND_", i, "_CELLCHAT"))
              )
            )
          )
        )
      )
    }
  )]
  LR_full[, paste0("RECEPTOR_", 1:3) := lapply(
    1:3,
    function(i) {
      ifelse(
        !is.na(get(paste0("RECEPTOR_", i, "_CELLPHONEDB"))),
        get(paste0("RECEPTOR_", i, "_CELLPHONEDB")),
        ifelse(
          !is.na(get(paste0("RECEPTOR_", i, "_ICELLNET"))),
          get(paste0("RECEPTOR_", i, "_ICELLNET")),
          ifelse(
            !is.na(get(paste0("RECEPTOR_", i, "_SCSR"))),
            get(paste0("RECEPTOR_", i, "_SCSR")),
            ifelse(
              !is.na(get(paste0("RECEPTOR_", i, "_SCTENSOR"))),
              get(paste0("RECEPTOR_", i, "_SCTENSOR")),
              ifelse(
                !is.na(get(paste0("RECEPTOR_", i, "_NICHENET"))),
                get(paste0("RECEPTOR_", i, "_NICHENET")),
                get(paste0("RECEPTOR_", i, "_CELLCHAT"))
              )
            )
          )
        )
      )
    }
  )]
  col_database <- paste0("DATABASE.1_", db_names)
  LR_full[, c("DATABASE") := do.call(paste, c(.SD, sep = ",")), .SDcols = col_database]
  LR_full[, c("DATABASE") := gsub("NA|NA,|,NA","", DATABASE)]
  LR_full[, c("N_DB") := rowSums(!is.na(.SD)), .SDcols = col_database]
  col_source <- paste0("SOURCE_", db_names)
  LR_full[, c("SOURCE") := do.call(paste, c(.SD, sep = ",")), .SDcols = col_source]
  LR_full[, c("SOURCE") := gsub("NA|NA,|,NA","", SOURCE)]
  col_anno <- paste0("ANNOTATION_", db_names)
  LR_full[, c("ANNOTATION") := do.call(paste, c(.SD, sep = ",")), .SDcols = col_anno]
  LR_full[, c("ANNOTATION") := gsub("NA|NA,|,NA","", ANNOTATION)]
  col_fam <- paste0("FAMILY_", db_names)
  LR_full[, c("FAMILY") := do.call(paste, c(.SD, sep = ",")), .SDcols = col_fam]
  LR_full[, c("FAMILY") := gsub("NA|NA,|,NA","", FAMILY)]
  col_subfam <- paste0("SUBFAMILY_", db_names)
  LR_full[, c("SUBFAMILY") := do.call(paste, c(.SD, sep = ",")), .SDcols = col_subfam]
  LR_full[, c("SUBFAMILY") := gsub("NA|NA,|,NA","", SUBFAMILY)]
  for(id_loop in c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3))) {
    cols_conf <- paste0(id_loop, "_CONF_", db_names)
    LR_full[, paste0(id_loop, "_CONF") := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = cols_conf]
    LR_full[, paste0(id_loop, "_CONF") := ifelse(is.na(get(paste0(id_loop, "_CONF"))) & !is.na(get(id_loop)),
                                                 1, get(paste0(id_loop, "_CONF")))]
    cols_type <- paste0(id_loop, "_TYPE_", db_names)
    LR_full[, paste0(id_loop, "_TYPE") := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = cols_type]
    LR_full[, paste0(id_loop, "_TYPE") := ifelse(is.na(get(paste0(id_loop, "_TYPE"))) & !is.na(get(id_loop)),
                                                 "provided", get(paste0(id_loop, "_TYPE")))]
  }
  rm_subLR <- sapply(LR_full$LR_SORTED, function(i) {
    sum(grepl(i, LR_full$LR_SORTED, fixed = TRUE))
  })
  LR_full[, keep_subLR := grepl("CELLPHONEDB|CELLCHAT|ICELLNET", DATABASE)]
  LR_full <- LR_full[rm_subLR == 1 | keep_subLR == TRUE]
  cols_to_keep <- c(
    paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3), "LR_SORTED",
    "DATABASE", "N_DB", "SOURCE", "ANNOTATION", "FAMILY", "SUBFAMILY",
    paste0("LIGAND_", 1:2, "_CONF"), paste0("LIGAND_", 1:2, "_TYPE"),
    paste0("RECEPTOR_", 1:3, "_CONF"), paste0("RECEPTOR_", 1:3, "_TYPE")
  )
  LR_full <- LR_full[, cols_to_keep, with = FALSE]


  if(curated) {
    LR_rm_sctensor <- c("SWISSPROT_STRING", "TREMBL_STRING")
    LR_rm_nichenet <- c("ppi_bidir_bidir", "ppi_bidir_bidir_go", "ppi_bidir_r",
                        "ppi_bidir_r_go", "ppi_l_bidir", "ppi_l_bidir_go",
                        "ppi_lr", "ppi_lr_go")
    LR_rm_scsr <- c("uniprot")
    LR_rm <- c(
      LR_rm_sctensor,
      LR_rm_nichenet,
      LR_rm_scsr,
      sapply(LR_rm_sctensor, function(i) {
        sapply(LR_rm_nichenet, function(j) {
          c(paste(i,j, sep = ","), paste(j, i, sep = ","))
        })
      }),
      sapply(LR_rm_sctensor, function(i) {
        sapply(LR_rm_scsr, function(j) {
          c(paste(i,j, sep = ","), paste(j, i, sep = ","))
        })
      }),
      sapply(LR_rm_scsr, function(i) {
        sapply(LR_rm_nichenet, function(j) {
          c(paste(i,j, sep = ","), paste(j, i, sep = ","))
        })
      }),
      sapply(LR_rm_sctensor, function(i) {
        sapply(LR_rm_nichenet, function(j) {
          sapply(LR_rm_scsr, function(k) {
            c(paste(i,j,k, sep = ","), paste(i,k,j, sep = ","), paste(j, i,k, sep = ","), paste(j, k, i, sep = ","),
              paste(k, i, j, sep = ","), paste(k, j, i, sep = ","))
          })
        })
      })
    )
    LR_full <- LR_full[!(SOURCE %in% LR_rm)]
  }
  return(LR_full)
}

#' Create a data.table with the ligand-receptor interactions from scTensor.
#'
#' @return data.table with ligand-receptor interactions and their relevant properties obtained from scTensor.
prepare_LR_scTensor <- function(
) {
  LR_SORTED <- LIGAND_1 <- RECEPTOR_1 <- NULL
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
    LR$LIGAND_1 <- LR_L_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  if(identical(LR_R_match$ENTREZID, as.character(LR$GENEID_R))) {
    LR$RECEPTOR_1 <- LR_R_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  data.table::setDT(LR)
  LR[, c("GENEID_L", "GENEID_R") := list(NULL, NULL)]
  data.table::setnames(LR, old = "SOURCEDB", new = "SOURCE")
  LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
    temp <- c(LIGAND_1[[i]], RECEPTOR_1[[i]])
    temp <- temp[!is.na(temp)]
    temp <- sort(temp)
    temp <- paste0(temp, collapse = "_")
  }))]
  LR <- LR[!duplicated(LR_SORTED)]
  cols_to_keep <- c(
    "LR_SORTED","SOURCE",
    "LIGAND_1", "RECEPTOR_1"
  )
  return(LR[, cols_to_keep, with = FALSE])
}

#' Create a data.table with the ligand-receptor interactions from SingleCellSignalR.
#'
#' @param one2one logical indicating if the orthology conversion has to be limited to one2one homology; default is FALSE.
#'
#' @return data.table with ligand-receptor interactions and their relevant properties obtained from SingleCellSignalR.
prepare_LR_scsr <- function(
  one2one = FALSE
) {
  LR <- SingleCellSignalR::LRdb
  data.table::setDT(LR)
  data.table::setnames(
    x = LR,
    old = c("ligand", "receptor"),
    new = c("L1", "R1")
  )
  ortho <- get_orthologs(
    genes = unique(c(LR$L1, LR$R1)),
    input_species = "human",
    one2one = one2one
  )
  ortho <- stats::na.omit(ortho)
  LR <- merge_LR_orthologs(
    LR_dt = LR,
    ortho_dt = ortho,
    nL = 1,
    charL = "L",
    nR = 1,
    charR = "R"
  )
  LR$SOURCE <-  ifelse(LR$PMIDs == "", LR$source, paste0(LR$source, ",PMID:", LR$PMIDs))

  cols_to_keep <- c(
    "LR_SORTED", "SOURCE",
    "LIGAND_1", "RECEPTOR_1",
    "LIGAND_1_CONF", "RECEPTOR_1_CONF",
    "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
  )
  return(LR[, cols_to_keep, with = FALSE])
}

#' #' Create a data.table with the ligand-receptor interactions from NICHENET.
#' #'
#' #' @param one2one logical indicating if the orthology conversion has to be limited to one2one homology; default is FALSE.
#' #'
#' #' @return data.table with ligand-receptor interactions and their relevant properties obtained from NICHENET.
#' prepare_LR_nichenet <- function(
#'   one2one = FALSE
#' ) {
#'   LR <- nichenetr::lr_network
#'   data.table::setDT(LR)
#'   data.table::setnames(
#'     x = LR,
#'     old = c("from", "to", "source"),
#'     new = c("L1", "R1", "SOURCE")
#'   )
#'   ortho <- get_orthologs(
#'     genes = unique(c(LR$L1, LR$R1)),
#'     input_species = "human",
#'     one2one = one2one
#'   )
#'   ortho <- stats::na.omit(ortho)
#'   LR <- merge_LR_orthologs(
#'     LR_dt = LR,
#'     ortho_dt = ortho,
#'     nL = 1,
#'     charL = "L",
#'     nR = 1,
#'     charR = "R"
#'   )
#'   cols_to_keep <- c(
#'     "LR_SORTED", "SOURCE",
#'     "LIGAND_1", "RECEPTOR_1",
#'     "LIGAND_1_CONF", "RECEPTOR_1_CONF",
#'     "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
#'   )
#'   return(LR[, cols_to_keep, with = FALSE])
#' }

#' Convert CellphoneDB database in a dataframe either complex or deconvoluted.
#'
#' @param deconvoluted logical indicating if returning the complex LR with 1 to 3 genes per category or the 1:1 deconvoluted data.frame.
#' @return data.table with ligand-receptor interactions and their relevant properties obtained from CELLPHONEDB.
create_LR_cpdb <- function(
  deconvoluted = FALSE
) {
  id_protein_multi_1_1 <- id_protein_multi_1_2 <- id_protein_1 <- id_protein_2 <- temp_id <-
     receptor_1 <- receptor_2 <- secreted_1 <- secreted_2 <- L_3 <- LR_ID <-  NULL
  data <- LRcp_raw
  dt_interac <- data.table::setDT(data$interaction_table)
  dt_multi <- data.table::setDT(data$multidata_table)
  dt_prot <- data.table::setDT(data$protein_table)
  dt_gene <- data.table::setDT(data$gene_table)
  dt_comp <- data.table::setDT(data$complex_composition_table)
  dt_full <- data.table::merge.data.table(
    x = dt_interac,
    y = dt_multi,
    by.x = "multidata_1_id",
    by.y = "id_multidata",
    all.x = TRUE,
    sort = FALSE
  )
  dt_full <- data.table::merge.data.table(
    x = dt_full,
    y = dt_multi,
    by.x = "multidata_2_id",
    by.y = "id_multidata",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_1", "_2")
  )
  dt_full <- data.table::merge.data.table(
    x = dt_full,
    y = dt_prot,
    by.x = "multidata_1_id",
    by.y = "protein_multidata_id",
    all.x = TRUE,
    sort = FALSE
  )
  dt_full <- data.table::merge.data.table(
    x = dt_full,
    y = dt_prot,
    by.x = "multidata_2_id",
    by.y = "protein_multidata_id",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_1", "_2")
  )
  dt_comp$id_comp <- paste0("id_protein_multi_", rowid(dt_comp$complex_multidata_id))
  dt_comp_dc <- dcast.data.table(
    dt_comp,
    formula = complex_multidata_id ~ id_comp,
    value.var = "protein_multidata_id"
  )
  dt_full <- data.table::merge.data.table(
    x = dt_full,
    y = dt_comp_dc,
    by.x = "multidata_1_id",
    by.y = "complex_multidata_id",
    all.x = TRUE,
    sort = FALSE
  )
  dt_full <- data.table::merge.data.table(
    x = dt_full,
    y = dt_comp_dc,
    by.x = "multidata_2_id",
    by.y = "complex_multidata_id",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_1", "_2")
  )
  dt_full[, id_protein_multi_1_1 := ifelse(is.na(id_protein_1), id_protein_multi_1_1, id_protein_1)]
  dt_full[, id_protein_multi_1_2 := ifelse(is.na(id_protein_2), id_protein_multi_1_2, id_protein_2)]
  dt_gene <- stats::na.omit(unique(dt_gene[, c("protein_id", "hgnc_symbol")]))
  dt_gene$temp_id <- data.table::rowid(dt_gene$protein_id)
  dt_gene <- dt_gene[temp_id == 1,]
  dt_gene[, temp_id := NULL]
  dt_full[, c(paste0("id_gene_multi_", 1:3, "_1"), paste0("id_gene_multi_", 1:3, "_2")) := c(
    lapply(1:3, function(i) {
      dt_gene[.SD, on = paste0("protein_id==id_protein_multi_", i, "_1"), get("x.hgnc_symbol")]
    }),
    lapply(1:3, function(i) {
      dt_gene[.SD, on = paste0("protein_id==id_protein_multi_", i, "_2"), get("x.hgnc_symbol")]
    })
  )
  ]
  dt_full[, paste0("L_", 1:3) := lapply(1:3, function(i) {
    ifelse(
      receptor_1 == 0 & receptor_2 == 1,
      get(paste0("id_gene_multi_", i, "_1")),
      ifelse(
        receptor_1 == 1 & receptor_2 == 0,
        get(paste0("id_gene_multi_", i, "_2")),
        ifelse(
          receptor_1 == 1 & receptor_2 == 1,
          ifelse(
            secreted_1 == 1 & secreted_2 == 0,
            get(paste0("id_gene_multi_", i, "_1")),
            ifelse(
              secreted_1 == 0 & secreted_2 == 1,
              get(paste0("id_gene_multi_", i, "_2")),
              get(paste0("id_gene_multi_", i, "_1"))
            )
          ),
          ifelse(
            secreted_1 == 1 & secreted_2 == 0,
            get(paste0("id_gene_multi_", i, "_1")),
            ifelse(
              secreted_1 == 0 & secreted_2 == 1,
              get(paste0("id_gene_multi_", i, "_2")),
              get(paste0("id_gene_multi_", i, "_1"))
            )
          )
        )
      )
    )
  })]
  dt_full[, paste0("R_", 1:3) := lapply(1:3, function(i) {
    ifelse(
      receptor_1 == 0 & receptor_2 == 1,
      get(paste0("id_gene_multi_", i, "_2")),
      ifelse(
        receptor_1 == 1 & receptor_2 == 0,
        get(paste0("id_gene_multi_", i, "_1")),
        ifelse(
          receptor_1 == 1 & receptor_2 == 1,
          ifelse(
            secreted_1 == 1 & secreted_2 == 0,
            get(paste0("id_gene_multi_", i, "_2")),
            ifelse(
              secreted_1 == 0 & secreted_2 == 1,
              get(paste0("id_gene_multi_", i, "_1")),
              get(paste0("id_gene_multi_", i, "_2"))
            )
          ),
          ifelse(
            secreted_1 == 1 & secreted_2 == 0,
            get(paste0("id_gene_multi_", i, "_2")),
            ifelse(
              secreted_1 == 0 & secreted_2 == 1,
              get(paste0("id_gene_multi_", i, "_1")),
              get(paste0("id_gene_multi_", i, "_2"))
            )
          )
        )
      )
    )
  })]
  dt_full <- dt_full[, c("id_cp_interaction", paste0("L_", 1:3), paste0("R_", 1:3))]
  dt_full[, L_3 := NULL]
  if(deconvoluted) {
    dt_full <- do.call("rbind", lapply(1:2, function(i) {
      do.call("rbind",lapply(1:3, function(j) {
        cols <- c(paste0("L_", i), paste0("R_",j))
        df <- stats::na.omit(dt_full[, cols, with = FALSE])
        colnames(df) <- c("L_1", "R_1")
        df$LR_ID <- paste(df$L_1, df$R_1, sep = "_")
        df <- df[!duplicated(df$LR_ID),]
      }) )
    }))
    dt_full <- dt_full[!duplicated(dt_full$LR_ID),]
    dt_full[, LR_ID := NULL]
  }
  return(dt_full)
}

#' Create a data.table with the ligand-receptor interactions from CELLPHONEDB.
#'
#' @param one2one logical indicating if the orthology conversion has to be limited to one2one homology; default is FALSE.
#' @param deconvoluted logical indiciating if the interactions are deconvoluted when complex; default is FALSE
#' @param keep_id logical indicating if keeping the column with the CELLPHONEDB interaction id
#'
#' @return data.table with ligand-receptor interactions and their relevant properties obtained from CELLPHONEDB.
prepare_LR_cpdb <- function(
  one2one = FALSE,
  deconvoluted = FALSE,
  keep_id = FALSE
) {
  LR <- create_LR_cpdb(
    deconvoluted = deconvoluted
  )
  genes_temp <- unique(unlist(LR))
  genes_temp <- genes_temp[!is.na(genes_temp)]
  ortho <- get_orthologs(
    genes = genes_temp,
    input_species = "human",
    one2one = one2one
  )
  ortho <- stats::na.omit(ortho)
  if(deconvoluted) {
    nL = 1
    nR = 1
    LR$SOURCE <-  "CPDB_DECONV"
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1",
      "LIGAND_1_CONF", "RECEPTOR_1_CONF",
      "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
    )
  } else {
    nL = 2
    nR = 3
    LR$SOURCE <-  "CPDB"
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
      "LIGAND_1_CONF", "LIGAND_2_CONF", "RECEPTOR_1_CONF", "RECEPTOR_2_CONF", "RECEPTOR_3_CONF",
      "LIGAND_1_TYPE", "LIGAND_2_TYPE", "RECEPTOR_1_TYPE", "RECEPTOR_2_TYPE", "RECEPTOR_3_TYPE"
    )
    if(keep_id){
      cols_to_keep <- c("id_cp_interaction", cols_to_keep)
    }
  }
  LR <- merge_LR_orthologs(
    LR_dt = LR,
    ortho_dt = ortho,
    nL = nL,
    charL = "L_",
    nR = nR,
    charR = "R_"
  )
  return(LR[, cols_to_keep, with = FALSE])
}

#' #' Create a data.table with the ligand-receptor interactions from CellChat.
#' #'
#' #' @return data.table with ligand-receptor interactions and their relevant properties obtained from CELLCHAT.
#' prepare_LR_CellChat <- function(
#' ) {
#'   LIGAND_1 <- RECEPTOR_1 <- RECEPTOR_2 <- LR_SORTED <- interaction_name_2 <- temp <- new <- NULL
#'   LR <- CellChat::CellChatDB.mouse$interaction
#'   setDT(LR)
#'   data.table::setnames(
#'     x = LR,
#'     old = c("evidence", "annotation"),
#'     new = c("SOURCE", "ANNOTATION")
#'   )
#'   LR[, LIGAND_1 := sub(" - .*", "", interaction_name_2) ]
#'   LR[, temp := sub(".* - ", "", interaction_name_2) ]
#'   LR[, RECEPTOR_1 := ifelse(grepl("+", temp, fixed = TRUE), gsub(".*\\((.+)\\+.*", "\\1", temp), temp)]
#'   LR[, RECEPTOR_2 := ifelse(grepl("+", temp, fixed = TRUE), gsub(".*\\+(.+)\\).*", "\\1", temp), NA)]
#'   LR[, temp := NULL]
#'   LR[, LIGAND_1 := gsub(" ", "", LIGAND_1)]
#'   LR[, RECEPTOR_1 := gsub(" ", "", RECEPTOR_1)]
#'   LR[, RECEPTOR_2 := gsub(" ", "", RECEPTOR_2)]
#'   #some CELLCHAT gene names (70) are not mgi_symbols and we need to convert them manually...
#'   convert_table <- CellChat_conversion
#'   genes_to_rm <- convert_table[new == "remove"]
#'   genes_to_change <- convert_table[new != "remove"]
#'   LR <- LR[!(LIGAND_1 %in% genes_to_rm$old) & !(RECEPTOR_1 %in% genes_to_rm$old) & !(RECEPTOR_2 %in% genes_to_rm$old)]
#'   LR[genes_to_change,
#'         `:=`(LIGAND_1 = new),
#'         on = "LIGAND_1==old"][
#'           genes_to_change,
#'           `:=`(RECEPTOR_1 = new),
#'           on = "RECEPTOR_1==old"][
#'             genes_to_change,
#'             `:=`(RECEPTOR_2 = new),
#'             on = "RECEPTOR_2==old"]
#'   LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
#'     temp <- c(LIGAND_1[[i]], RECEPTOR_1[[i]], RECEPTOR_2[[i]])
#'     temp <- temp[!is.na(temp)]
#'     temp <- sort(temp)
#'     temp <- paste0(temp, collapse = "_")
#'   }))]
#'   LR <- LR[!duplicated(LR_SORTED)]
#'   cols_to_keep <- c(
#'     "LR_SORTED",
#'     "ANNOTATION", "SOURCE",
#'     "LIGAND_1", "RECEPTOR_1", "RECEPTOR_2"
#'   )
#'   LR <- LR[, cols_to_keep, with = FALSE]
#'   return(LR)
#' }

#' Create a data.table with the ligand-receptor interactions from ICELLNET.
#'
#' @param one2one logical indicating if the orthology conversion has to be limited to one2one homology; default is FALSE.
#'
#' @return data.table with ligand-receptor interactions and their relevant properties obtained from ICELLNET.
prepare_LR_ICELLNET <- function(
  one2one = FALSE
) {
  LR <- utils::read.csv(
    file = curl::curl(url = "https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/database.tsv"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = ""
  )
  data.table::setDT(LR)
  data.table::setnames(
    x = LR,
    old = c("Ligand 1", "Ligand 2",
            "Receptor 1", "Receptor 2", "Receptor 3",
            "PubMed ID", "Family", "Subfamily", "Classifications"),
    new = c("L1", "L2",
            "R1", "R2", "R3",
            "SOURCE", "FAMILY", "SUBFAMILY", "ANNOTATION")
  )
  ortho <- get_orthologs(
    genes = unique(c(LR$L1, LR$L2, LR$R1, LR$R2, LR$R3)),
    input_species = "human",
    one2one = one2one
  )
  ortho <- stats::na.omit(ortho)
  LR <- merge_LR_orthologs(
    LR_dt = LR,
    ortho_dt = ortho,
    nL = 2,
    charL = "L",
    nR = 3,
    charR = "R"
  )
  cols_to_keep <- c(
    "LR_SORTED",
    "FAMILY", "SUBFAMILY", "ANNOTATION", "SOURCE",
    "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
    "LIGAND_1_CONF", "LIGAND_2_CONF", "RECEPTOR_1_CONF", "RECEPTOR_2_CONF", "RECEPTOR_3_CONF",
    "LIGAND_1_TYPE", "LIGAND_2_TYPE", "RECEPTOR_1_TYPE", "RECEPTOR_2_TYPE", "RECEPTOR_3_TYPE"
  )
  return(LR[, cols_to_keep, with = FALSE])
}

#' Take a LR data.table of human genes and convert it with mouse genes.
#'
#' @param LR_dt data.table of human ligand-recetpor interactions.
#' @param ortho_dt data.table of human-mouse homologues.
#' @param nL numeric indicating the number of ligand columns in LR_dt.
#' @param charL character indicating the prefix name of the ligand columns in LR_dt.
#' @param nR numeric indicating the number of receptor columns in LR_dt.
#' @param charR character indicating the prefix name of the receptor columns in LR_dt.
#'
#' @return data.table with the mouse ligand-receptor interactions. Filtering of duplicated is also performed.
merge_LR_orthologs <- function(
  LR_dt,
  ortho_dt,
  nL,
  charL,
  nR,
  charR
) {
  to_keep <- LR_SORTED <- NULL
  out_names <- c(
    sapply(1:nL, function(i) {paste0("LIGAND_", i, c("", "_CONF", "_TYPE"))}),
    sapply(1:nR, function(i) {paste0("RECEPTOR_", i, c("", "_CONF", "_TYPE"))})
  )
  merge_id <- c("mouse_symbol", "confidence", "type")
  LR_temp <- data.table::copy(LR_dt)
  LR_temp[, c(out_names) :=
       c(
         sapply(
           1:nL,
           function(i) {
             as.list(
               ortho_dt[.SD,
                     on = c(paste0("human_symbol==", charL, i)),
                     mget(paste0("x.", merge_id))
                     ])
           }
         ),
         sapply(
           1:nR,
           function(i) {
             as.list(
               ortho_dt[.SD,
                     on = c(paste0("human_symbol==", charR, i)),
                     mget(paste0("x.", merge_id))
                     ])
           }
         )
       )
     ]
  LR_temp <- stats::na.omit(LR_temp, cols = c("LIGAND_1", "RECEPTOR_1"))
  LR_temp[, to_keep := sapply(1:nrow(.SD), function(i) {
    all(c(sapply(1:nL, function(j) {!(is.na(get(paste0("LIGAND_", j))[[i]]) & !is.na(get(paste0(charL, j))[[i]]))}),
          sapply(1:nR, function(j) {!(is.na(get(paste0("RECEPTOR_", j))[[i]]) & !is.na(get(paste0(charR, j))[[i]]))})))
  })]
  LR_temp <- LR_temp[to_keep == TRUE]
  LR_temp[, to_keep := NULL]
  LR_temp[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
    temp <- c(sapply(1:nL, function(j) {get(paste0("LIGAND_", j))[[i]]}),
              sapply(1:nR, function(j) {get(paste0("RECEPTOR_", j))[[i]]}))
    temp <- temp[!is.na(temp)]
    temp <- sort(temp)
    temp <- paste0(temp, collapse = "_")
  }))]
  LR_temp <- LR_temp[!duplicated(LR_SORTED)]
  return(LR_temp)
}

get_orthologs <- function(
  genes,
  input_species,
  one2one = FALSE
) {
  ensembl_gene_id <- inl <- outl <- output <- input <- confidence <- NULL
  if(input_species == "mouse") {
    id_in <- "mmusculus"
    id_out <- "hsapiens"
    id_gene <- "mgi_symbol"
    name_in <- "mouse"
    name_out <- "human"
  } else if(input_species == "human") {
    id_in <- "hsapiens"
    id_out <- "mmusculus"
    id_gene <- "hgnc_symbol"
    name_in <- "human"
    name_out <- "mouse"
  } else {
    stop("Species not supported in function get_orthologs.")
  }
  dataset <- paste0(id_in, "_gene_ensembl")
  gene_name <- paste0(id_out, "_homolog_associated_gene_name")
  ortho_confidence <- paste0(id_out, "_homolog_orthology_confidence")
  ortho_type <- paste0(id_out, "_homolog_orthology_type")
  mart <- biomaRt::useMart(
    "ensembl",
    dataset = dataset
  )
  ensembl <- biomaRt::getBM(
    attributes = c(
      id_gene,
      "ensembl_gene_id"
    ),
    filters = id_gene,
    mart = mart,
    values = genes
  )
  ensembl_conv <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      gene_name,
      ortho_confidence,
      ortho_type
    ),
    filters = "ensembl_gene_id",
    mart = mart,
    value = ensembl$ensembl_gene_id)
  data.table::setDT(ensembl)
  data.table::setDT(ensembl_conv)
  ensembl_all <- data.table::merge.data.table(
    x = ensembl,
    y = ensembl_conv,
    by = "ensembl_gene_id",
    all = TRUE,
    sort = FALSE
  )
  if(one2one) {
    ensembl_all <- ensembl_all[eval(as.symbol(ortho_type)) == 'ortholog_one2one',]
  } else {
    #ensembl_all <- ensembl_all[eval(as.symbol(ortho_type)) %in% c('ortholog_one2one','ortholog_one2many'),]
  }
  ensembl_all <- ensembl_all[, ensembl_gene_id := NULL]
  ensembl_all <- unique(ensembl_all)
  data.table::setnames(
    x = ensembl_all,
    old = c(id_gene, gene_name, ortho_confidence, ortho_type),
    new = c("input", "output", "confidence", "type")
  )
  ensembl_all <- stats::na.omit(ensembl_all)
  #check for remaining duplicate
  if(sum(duplicated(ensembl_all[["input"]])) > 0) {
    ensembl_all[, inl := tolower(input)]
    ensembl_all[, outl := tolower(output)]
    dup_input <- unique(ensembl_all$input[duplicated(ensembl_all$input)])
    for(g in dup_input) {
      dt_g <- ensembl_all[input == g]
      dt_gconf <- dt_g[confidence == 1]
      if(nrow(dt_gconf) == 1) {
        ensembl_all <- ensembl_all[!(input == g & confidence == 0)]
      } else if(nrow(dt_gconf) == 0) {
        dt_gsame <- dt_g[inl == outl]
        if(nrow(dt_gsame) == 0) {
          g_keep <- dt_g[1]$output
        } else {
          g_keep <- dt_gsame[1]$output
        }
        ensembl_all <- ensembl_all[!(input == g & output != g_keep)]
      } else {
        dt_gsame <- dt_gconf[inl == outl]
        if(nrow(dt_gsame) == 0) {
          g_keep <- dt_gconf[1]$output
        } else {
          g_keep <- dt_gsame[1]$output
        }
        ensembl_all <- ensembl_all[!(input == g & confidence == 0)]
        ensembl_all <- ensembl_all[!(input == g & output != g_keep)]
      }
    }
    ensembl_all[, inl := NULL]
    ensembl_all[, outl := NULL]
    if(sum(duplicated(ensembl_all[["input"]])) > 0) {
      warning("There are some duplicates from orthology conversion. Removing them by using 'unique'.")
      ensembl_all <- unique(ensembl_all, by = "input")
    }
  }
  data.table::setnames(
    x = ensembl_all,
    old = c("input", "output"),
    new = c(paste0(name_in, "_symbol"), paste0(name_out, "_symbol"))
  )
  return(ensembl_all)
}

