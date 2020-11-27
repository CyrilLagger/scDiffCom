#' Build a data.table of curated ligand-receptor interactions obtained from 6 databases.
#'
#' @param species human or mouse
#'
#' @return A data.table with ligands, receptors and some annotations (database of origin and source of curation).
build_LRdb <- function(
  species
) {
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  LRdb_notCurated <- combine_LR_db(
    species = species,
    one2one = FALSE,
    curated = FALSE
  )
  LRdb_curated <- combine_LR_db(
    species = species,
    one2one = FALSE,
    curated = TRUE
  )
  LRdb_GO <- get_GO_interactions(
    LR_db = LRdb_curated
  )
  return(list(
    LRdb_notCurated = LRdb_notCurated,
    LRdb_curated = LRdb_curated,
    LRdb_curated_GO = LRdb_GO
  ))
}

get_GO_interactions <- function(
  LR_db
) {
  GO_NAME <- NULL
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package \"biomaRt\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("ontoProc", quietly = TRUE)) {
    stop("Package \"ontoProc\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("ontologyIndex", quietly = TRUE)) {
    stop("Package \"ontologyIndex\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  mgi_symbol <- NULL
  LR_genes <- unique(unlist(LR_db[, c("LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3")]))
  LR_genes <- LR_genes[!is.na(LR_genes)]
  mart <- biomaRt::useMart(
    "ensembl",
    dataset = "mmusculus_gene_ensembl"
  )
  LR_genes_info <- biomaRt::getBM(
    attributes = c(
      "mgi_symbol",
      "go_id",
      "name_1006"
    ),
    filters = "mgi_symbol",
    mart = mart,
    values = LR_genes
  )
  setDT(LR_genes_info)
  onto_go_terms <- ontoProc::getGeneOnto()
  go_names <- onto_go_terms$name
  LR_genes_go <- sapply(
    LR_genes,
    function(gene) {
      temp_go <- unique(LR_genes_info[mgi_symbol == gene]$name_1006)
      ontologyIndex::get_ancestors(onto_go_terms, names(go_names[go_names %in% temp_go]))
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
  # LR_interactions_go_union <- rbindlist(
  #   apply(
  #     LR_db,
  #     MARGIN = 1,
  #     function(row) {
  #       LIGAND_GO <- unique(c(
  #         LR_genes_go[[row[["LIGAND_1"]]]],
  #         LR_genes_go[[row[["LIGAND_2"]]]]
  #       ))
  #       RECEPTOR_GO <- unique(c(
  #         LR_genes_go[[row[["RECEPTOR_1"]]]],
  #         LR_genes_go[[row[["RECEPTOR_2"]]]],
  #         LR_genes_go[[row[["RECEPTOR_3"]]]]
  #       ))
  #       res_union <- unique(c(LIGAND_GO, RECEPTOR_GO))
  #       res_union <- data.table(
  #         LR_SORTED = rep(row[["LR_SORTED"]], length(res_union)),
  #         GO_ID = res_union
  #       )
  #     }
  #   )
  # )
  LR_interactions_go_intersection <- rbindlist(
    apply(
      LR_db,
      MARGIN = 1,
      function(row) {
        LIGAND_GO <- unique(c(
          LR_genes_go[[row[["LIGAND_1"]]]],
          LR_genes_go[[row[["LIGAND_2"]]]]
        ))
        RECEPTOR_GO <- unique(c(
          LR_genes_go[[row[["RECEPTOR_1"]]]],
          LR_genes_go[[row[["RECEPTOR_2"]]]],
          LR_genes_go[[row[["RECEPTOR_3"]]]]
        ))
        res_inter <- intersect(LIGAND_GO, RECEPTOR_GO)
        if (length(res_inter) > 0) {
          res_inter <- data.table(
            LR_SORTED = rep(row[["LR_SORTED"]], length(res_inter)),
            GO_ID = res_inter
          )
        } else {
          res_inter <- NULL
        }
        return(res_inter)
      }
    )
  )
  go_id_name_dt <- data.table(
    GO_name = go_names,
    ID = names(go_names)
  )
  # LR_interactions_go_union[
  #   go_id_name_dt,
  #   on = "GO_ID==ID",
  #   GO_NAME := i.GO_name
  #   ]
  LR_interactions_go_intersection[
    go_id_name_dt,
    on = "GO_ID==ID",
    GO_NAME := i.GO_name
    ]
  return(list(
    # LR_GO_union = LR_interactions_go_union,
    LR_GO_intersection = LR_interactions_go_intersection
  ))
}

combine_LR_db <- function(
  species,
  one2one = FALSE,
  curated = TRUE
) {
  DATABASE <- SOURCE <- ANNOTATION <- FAMILY <- SUBFAMILY <- keep_subLR <-
    SOURCE_CLEAN <- SOURCE_no_digit <- is_complex_temp <- LIGAND_2 <- RECEPTOR_2 <- LR_vectorized_temp <-
    N_IS_SUBPART <- NULL
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  #already curated
  LR_CONNECTOMEDB <- prepare_LR_CONNECTOMEDB(species = species, one2one = one2one) #PMID:
  LR_celltalk <- prepare_LR_CellTalkDB(species = species) #PMID:
  LR_ic <- prepare_LR_ICELLNET(species = species, one2one = one2one) #PMID:
  LR_cpdb <- prepare_LR_cpdb(species = species, one2one = one2one, deconvoluted = FALSE) #CPDB
  LR_cc <- prepare_LR_CellChat(species = species) #PMID: KEGG: PMC: PMC
  #not fully curated
  LR_scsr <- prepare_LR_scsr(species = species, one2one = one2one) #"PMID:" "Ramilowski2015" "HPMR" "HPRD" and
  #"reactome" "IUPHAR" "PPI" "cellsignal.com"
  LR_niche <- prepare_LR_nichenet(species = species, one2one = one2one) #"KEGG:nichenet"  "IUPHAR" "Ramilowski2015" "PPI"
  LR_sct <- prepare_LR_scTensor(species = species) #"SCT:DLRP" "IUPHAR" "HPMR" "SCT:CPDB" "SCT:SCSR"   "PPI"
  #"HPRD" SCT:Ramilowski2015"
  if (curated) {
    LR_scsr <- LR_scsr[SOURCE != "PPI"]
    LR_niche <- LR_niche[SOURCE != "PPI"]
    LR_sct <- LR_sct[SOURCE != "PPI"]
  }
  LR_full <- rbindlist(
    list(
      "CONNECTOMEDB" = LR_CONNECTOMEDB,
      "CELLTALK" = LR_celltalk,
      "SCTENSOR" = LR_sct,
      "SCSR" = LR_scsr,
      "NICHENET" = LR_niche,
      "CELLPHONEDB" = LR_cpdb,
      "CELLCHAT" = LR_cc,
      "ICELLNET" = LR_ic
    ),
    use.names = TRUE,
    fill = TRUE,
    idcol = "DATABASE"
  )
  col_md <- colnames(LR_full)
  col_md <- col_md[col_md != "LR_SORTED"]
  db_names <- c("CONNECTOMEDB", "CELLTALK", "CELLCHAT", "CELLPHONEDB", "SCSR", "SCTENSOR", "ICELLNET", "NICHENET")
  LR_full <- dcast.data.table(
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
                ifelse(
                  !is.na(get(paste0("LIGAND_", i, "_CONNECTOMEDB"))),
                  get(paste0("LIGAND_", i, "_CONNECTOMEDB")),
                  ifelse(
                    !is.na(get(paste0("LIGAND_", i, "_CELLTALK"))),
                    get(paste0("LIGAND_", i, "_CELLTALK")),
                    get(paste0("LIGAND_", i, "_CELLCHAT"))
                  )
                )
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
                ifelse(
                  !is.na(get(paste0("RECEPTOR_", i, "_CONNECTOMEDB"))),
                  get(paste0("RECEPTOR_", i, "_CONNECTOMEDB")),
                  ifelse(
                    !is.na(get(paste0("RECEPTOR_", i, "_CELLTALK"))),
                    get(paste0("RECEPTOR_", i, "_CELLTALK")),
                    get(paste0("RECEPTOR_", i, "_CELLCHAT"))
                  )
                )
              )
            )
          )
        )
      )
    }
  )]
  col_database <- paste0("DATABASE.1_", db_names)
  LR_full[, c("DATABASE") := do.call(paste, c(.SD, sep = ";")), .SDcols = col_database]
  LR_full[, c("DATABASE") := gsub("NA|NA;|;NA", "", DATABASE)]
  LR_full[, c("N_DB") := rowSums(!is.na(.SD)), .SDcols = col_database]
  col_source <- paste0("SOURCE_", db_names)
  LR_full[, c("SOURCE") := do.call(paste, c(.SD, sep = ";")), .SDcols = col_source]
  LR_full[, c("SOURCE") := gsub("NA|NA;|;NA", "", SOURCE)]
  col_anno <- paste0("ANNOTATION_", db_names)
  LR_full[, c("ANNOTATION") := do.call(paste, c(.SD, sep = ";")), .SDcols = col_anno]
  LR_full[, c("ANNOTATION") := gsub("NA|NA;|;NA", "", ANNOTATION)]
  col_fam <- paste0("FAMILY_", db_names)
  LR_full[, c("FAMILY") := do.call(paste, c(.SD, sep = ";")), .SDcols = col_fam]
  LR_full[, c("FAMILY") := gsub("NA|NA;|;NA", "", FAMILY)]
  col_subfam <- paste0("SUBFAMILY_", db_names)
  LR_full[, c("SUBFAMILY") := do.call(paste, c(.SD, sep = ";")), .SDcols = col_subfam]
  LR_full[, c("SUBFAMILY") := gsub("NA|NA;|;NA", "", SUBFAMILY)]
  if(species == "mouse") {
    for (id_loop in c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3))) {
      cols_conf <- paste0(id_loop, "_CONF_", db_names)
      LR_full[, paste0(id_loop, "_CONF") := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = cols_conf]
      LR_full[, paste0(id_loop, "_CONF") := ifelse(is.na(get(paste0(id_loop, "_CONF"))) & !is.na(get(id_loop)),
                                                   1, get(paste0(id_loop, "_CONF"))
      )]
      cols_type <- paste0(id_loop, "_TYPE_", db_names)
      LR_full[, paste0(id_loop, "_TYPE") := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = cols_type]
      LR_full[, paste0(id_loop, "_TYPE") := ifelse(is.na(get(paste0(id_loop, "_TYPE"))) & !is.na(get(id_loop)),
                                                   "provided", get(paste0(id_loop, "_TYPE"))
      )]
    }
  }
  if (curated) {
    LR_full[, is_complex_temp := fifelse(!is.na(LIGAND_2) | !is.na(RECEPTOR_2), TRUE, FALSE)]
    LR_full[, LR_vectorized_temp := list(sapply(1:nrow(.SD), function(i) {
      temp <- c(
        sapply(1:2, function(j) {
          get(paste0("LIGAND_", j))[[i]]
        }),
        sapply(1:3, function(j) {
          get(paste0("RECEPTOR_", j))[[i]]
        })
      )
      temp <- temp[!is.na(temp)]
    }))]
    LR_full_simple <- LR_full[is_complex_temp == FALSE]
    LR_full_complex <- LR_full[is_complex_temp == TRUE]
    LR_full_simple_list <- LR_full_simple$LR_vectorized_temp
    LR_full_complex_list <- LR_full_complex$LR_vectorized_temp
    rm_LR_full <- sapply(
      LR_full_simple_list,
      function(i) {
        sum(
          sapply(
            LR_full_complex_list,
            function(j) {
              all(i %in% j)
            }
          )
        )
      }
    )
    LR_full_simple[, N_IS_SUBPART := rm_LR_full]
    LR_full[LR_full_simple, on = "LR_SORTED", N_IS_SUBPART := i.N_IS_SUBPART]
    setnafill(LR_full, fill = 0, cols = "N_IS_SUBPART")
    LR_full[, keep_subLR := grepl("CELLPHONEDB|CELLCHAT|ICELLNET", DATABASE)]
    LR_full <- LR_full[N_IS_SUBPART == 0 | keep_subLR == TRUE]
    LR_full <- LR_full[DATABASE != "SCTENSOR"]
    cols_to_keep <- NULL
  } else {
    cols_to_keep <- NULL
  }
  if (species == "mouse") {
    cols_to_keep <- c(
      cols_to_keep,
      paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3), "LR_SORTED",
      "DATABASE", "N_DB", "SOURCE", "ANNOTATION", "FAMILY", "SUBFAMILY",
      paste0("LIGAND_", 1:2, "_CONF"), paste0("LIGAND_", 1:2, "_TYPE"),
      paste0("RECEPTOR_", 1:3, "_CONF"), paste0("RECEPTOR_", 1:3, "_TYPE")
    )
  }
  if (species == "human") {
    cols_to_keep <- c(
      cols_to_keep,
      paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3), "LR_SORTED",
      "DATABASE", "N_DB", "SOURCE", "ANNOTATION", "FAMILY", "SUBFAMILY"
    )
  }
  LR_full <- LR_full[, cols_to_keep, with = FALSE]
  return(LR_full)
}

prepare_LR_CONNECTOMEDB <- function(
  species,
  one2one = FALSE
) {
  LR_SORTED <- SOURCE <- NULL
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  LR <- NATMI_connectomeDB2020
  setDT(LR)
  setnames(
    x = LR,
    old = c("Ligand.gene.symbol", "Receptor.gene.symbol", "PMID.support"),
    new = c("L1", "R1", "SOURCE")
  )
  if (species == "mouse") {
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
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1",
      "LIGAND_1_CONF", "RECEPTOR_1_CONF",
      "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
    )
  }
  if (species == "human") {
    setnames(
      x = LR,
      old = c(
        "L1", "R1"
      ),
      new = c(
        "LIGAND_1", "RECEPTOR_1"
      )
    )
    LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
      temp <- c(
        sapply(1:1, function(j) {
          get(paste0("LIGAND_", j))[[i]]
        }),
        sapply(1:1, function(j) {
          get(paste0("RECEPTOR_", j))[[i]]
        })
      )
      temp <- temp[!is.na(temp)]
      temp <- sort(temp)
      temp <- paste0(temp, collapse = "_")
    }))]
    LR <- LR[!duplicated(LR_SORTED)]
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1"
    )
  }
  LR[, SOURCE := gsub("\\s", "", SOURCE)]
  LR[, SOURCE := paste0("PMID:", SOURCE)]
  LR[, SOURCE := gsub(",", ";PMID:", SOURCE)]
  return(LR[, cols_to_keep, with = FALSE])
}

prepare_LR_CellTalkDB <- function(
  species
) {
  LR_SORTED <- SOURCE <- NULL
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  if (species == "mouse") {
    LR <- CellTalkDB_mouse
  }
  if (species == "human") {
    LR <- CellTalkDB_human
  }
  setDT(LR)
  setnames(
    LR,
    old = c("ligand_gene_symbol", "receptor_gene_symbol", "evidence"),
    new = c("LIGAND_1", "RECEPTOR_1", "SOURCE")
  )
  LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
    temp <- c(
      sapply(1:1, function(j) {
        get(paste0("LIGAND_", j))[[i]]
      }),
      sapply(1:1, function(j) {
        get(paste0("RECEPTOR_", j))[[i]]
      })
    )
    temp <- temp[!is.na(temp)]
    temp <- sort(temp)
    temp <- paste0(temp, collapse = "_")
  }))]
  LR <- LR[!duplicated(LR_SORTED)]
  cols_to_keep <- c(
    "LR_SORTED", "SOURCE",
    "LIGAND_1", "RECEPTOR_1"
  )
  LR[, SOURCE := paste0("PMID:", SOURCE)]
  LR[, SOURCE := gsub(",", ";PMID:", SOURCE)]
  return(LR[, cols_to_keep, with = FALSE])
}

prepare_LR_scTensor <- function(
  species
) {
  SOURCE <- NULL
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package \"AnnotationDbi\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("LRBase.Mmu.eg.db", quietly = TRUE)) {
    stop("Package \"LRBase.Mmu.eg.db\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("LRBase.Hsa.eg.db", quietly = TRUE)) {
    stop("Package \"LRBase.Hsa.eg.db\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
    stop("Package \"AnnotationHub\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  LR_SORTED <- LIGAND_1 <- RECEPTOR_1 <- NULL
  if (species == "mouse") {
    key <- AnnotationDbi::keys(
      LRBase.Mmu.eg.db::LRBase.Mmu.eg.db,
      keytype = "GENEID_L"
    )
    LR <- AnnotationDbi::select(
      LRBase.Mmu.eg.db::LRBase.Mmu.eg.db,
      keys = key,
      columns = c("GENEID_L", "GENEID_R", "SOURCEDB"),
      keytype = "GENEID_L"
    )
  }
  if (species == "human"){
    key <- AnnotationDbi::keys(
      LRBase.Hsa.eg.db::LRBase.Hsa.eg.db,
      keytype = "GENEID_L"
    )
    LR <- AnnotationDbi::select(
      LRBase.Hsa.eg.db::LRBase.Hsa.eg.db,
      keys = key,
      columns = c("GENEID_L", "GENEID_R", "SOURCEDB"),
      keytype = "GENEID_L"
    )
  }
  LR <- unique(LR)
  ah <- AnnotationHub::AnnotationHub()
  if (species == "mouse") {
    hs <- AnnotationHub::query(
      ah,
      c("OrgDb", "Mus musculus")
    )[[1]]
  }
  if (species == "human") {
    hs <- AnnotationHub::query(
      ah,
      c("OrgDb", "Homo sapiens")
    )[[1]]
  }
  LR_L_match <- AnnotationDbi::select(
    hs,
    column = c("SYMBOL", "ENTREZID"),
    keytype = "ENTREZID",
    keys = as.character(LR$GENEID_L)
  )
  LR_R_match <- AnnotationDbi::select(
    hs,
    column = c("SYMBOL", "ENTREZID"),
    keytype = "ENTREZID",
    keys = as.character(LR$GENEID_R)
  )
  if (identical(LR_L_match$ENTREZID, as.character(LR$GENEID_L))) {
    LR$LIGAND_1 <- LR_L_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  if (identical(LR_R_match$ENTREZID, as.character(LR$GENEID_R))) {
    LR$RECEPTOR_1 <- LR_R_match$SYMBOL
  } else {
    stop("Matching not possible.")
  }
  setDT(LR)
  LR[, c("GENEID_L", "GENEID_R") := list(NULL, NULL)]
  setnames(LR, old = "SOURCEDB", new = "SOURCE")
  LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
    temp <- c(LIGAND_1[[i]], RECEPTOR_1[[i]])
    temp <- temp[!is.na(temp)]
    temp <- sort(temp)
    temp <- paste0(temp, collapse = "_")
  }))]
  LR <- LR[!duplicated(LR_SORTED)]
  LR <- LR[!is.na(LIGAND_1) & !is.na(RECEPTOR_1)]
  cols_to_keep <- c(
    "LR_SORTED", "SOURCE",
    "LIGAND_1", "RECEPTOR_1"
  )
  if (species == "mouse") {
    LR[, SOURCE := gsub("ENSEMBL_DLRP|NCBI_DLRP", "SCT:DLRP", SOURCE)]
    LR[, SOURCE := gsub("ENSEMBL_IUPHAR|NCBI_IUPHAR", "IUPHAR", SOURCE)]
    LR[, SOURCE := gsub("ENSEMBL_HPMR|NCBI_HPMR", "HPMR", SOURCE)]
    LR[, SOURCE := gsub("ENSEMBL_CELLPHONEDB|NCBI_CELLPHONEDB", "SCT:CPDB", SOURCE)]
    LR[, SOURCE := gsub("ENSEMBL_SINGLECELLSIGNALR|NCBI_SINGLECELLSIGNALR", "SCT:SCSR", SOURCE)]
    LR[, SOURCE := gsub("SWISSPROT_STRING", "PPI", SOURCE)]
    LR[, SOURCE := gsub("TREMBL_STRING", "PPI", SOURCE)]
  }
  if (species == "human") {
    LR[, SOURCE := gsub("SWISSPROT_STRING|TREMBL_STRING|BADERLAB", "PPI", SOURCE)]
    LR[, SOURCE := gsub("SWISSPROT_HPRD|TREMBL_HPRD", "HPRD", SOURCE)]
    LR[, SOURCE := gsub("IUPHAR", "IUPHAR", SOURCE)]
    LR[, SOURCE := gsub("DLRP", "SCT:DLRP", SOURCE)]
    LR[, SOURCE := gsub("HPMR", "HPMR", SOURCE)]
    LR[, SOURCE := gsub("FANTOM5", "SCT:Ramilowski2015", SOURCE)]
    LR[, SOURCE := gsub("CELLPHONEDB", "SCT:CPDB", SOURCE)]
    LR[, SOURCE := gsub("SINGLECELLSIGNALR", "SCT:SCSR", SOURCE)]
  }
  return(LR[, cols_to_keep, with = FALSE])
}

prepare_LR_scsr <- function(
  species,
  one2one = FALSE
) {
  LR_SORTED <- SOURCE <- PMIDs <- NULL
  if (!requireNamespace("SingleCellSignalR", quietly = TRUE)) {
    stop("Package \"SingleCellSignalR\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  LR <- SingleCellSignalR::LRdb
  setDT(LR)
  setnames(
    x = LR,
    old = c("ligand", "receptor"),
    new = c("L1", "R1")
  )
  if (species == "mouse") {
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
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1",
      "LIGAND_1_CONF", "RECEPTOR_1_CONF",
      "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
    )
  }
  if (species == "human") {
    setnames(
      x = LR,
      old = c(
        "L1", "R1"
      ),
      new = c(
        "LIGAND_1", "RECEPTOR_1"
      )
    )
    LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
      temp <- c(
        sapply(1:1, function(j) {
          get(paste0("LIGAND_", j))[[i]]
        }),
        sapply(1:1, function(j) {
          get(paste0("RECEPTOR_", j))[[i]]
        })
      )
      temp <- temp[!is.na(temp)]
      temp <- sort(temp)
      temp <- paste0(temp, collapse = "_")
    }))]
    LR <- LR[!duplicated(LR_SORTED)]
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1"
    )
  }
  LR[, PMIDs := ifelse(PMIDs == "", NA, PMIDs )]
  LR[, PMIDs := ifelse(is.na(PMIDs), NA, paste0("PMID:", PMIDs))]
  LR[, PMIDs := gsub(",", ";PMID:", PMIDs)]
  LR[, PMIDs := gsub(")", "", PMIDs)]
  LR[, PMIDs := gsub(" ", "", PMIDs)]
  LR[, source := gsub(",", ";", source)]
  LR[, source := gsub("fantom5", "Ramilowski2015", source)]
  LR[, source := gsub("literature;", "", source)]
  LR[, source := gsub("literature", "", source)]
  LR[, source := gsub("uniprot", "PPI", source)]
  LR[, source := ifelse(source == "", NA, source)]
  LR[, SOURCE := paste(PMIDs, source, sep = ";")]
  LR[, SOURCE := gsub(";NA", "", SOURCE)]
  LR[, SOURCE := gsub("NA;", "", SOURCE)]
  return(LR[, cols_to_keep, with = FALSE])
}

prepare_LR_nichenet <- function(
  species,
  one2one = FALSE
) {
  LR_SORTED <- SOURCE <- R3 <- NULL
  if (!requireNamespace("nichenetr", quietly = TRUE)) {
    stop("Package \"nichenetr\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  LR <- nichenetr::lr_network
  setDT(LR)
  setnames(
    x = LR,
    old = c("from", "to", "source"),
    new = c("L1", "R1", "SOURCE")
  )
  if (species == "mouse") {
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
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1",
      "LIGAND_1_CONF", "RECEPTOR_1_CONF",
      "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
    )
  }
  if (species == "human") {
    setnames(
      x = LR,
      old = c(
        "L1", "R1"
      ),
      new = c(
        "LIGAND_1", "RECEPTOR_1"
      )
    )
    LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
      temp <- c(
        sapply(1:1, function(j) {
          get(paste0("LIGAND_", j))[[i]]
        }),
        sapply(1:1, function(j) {
          get(paste0("RECEPTOR_", j))[[i]]
        })
      )
      temp <- temp[!is.na(temp)]
      temp <- sort(temp)
      temp <- paste0(temp, collapse = "_")
    }))]
    LR <- LR[!duplicated(LR_SORTED)]
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1"
    )
  }
  source_temp <- unique(LR$SOURCE)
  source_temp_kegg <- paste0(source_temp[grepl("kegg", source_temp)], collapse = "|")
  source_temp_ppi <- paste0(source_temp[grepl("ppi", source_temp)], collapse = "|")
  LR[, SOURCE := gsub(source_temp_ppi, "PPI", SOURCE)]
  LR[, SOURCE := gsub(source_temp_kegg, "KEGG:nichenet", SOURCE)]
  LR[, SOURCE := gsub("pharmacology", "IUPHAR", SOURCE)]
  LR[, SOURCE := gsub("ramilowski_known", "Ramilowski2015", SOURCE)]
  return(LR[, cols_to_keep, with = FALSE])
}

create_LR_cpdb <- function(
  deconvoluted = FALSE
) {
  id_protein_multi_1_1 <- id_protein_multi_1_2 <- id_protein_1 <- id_protein_2 <- temp_id <-
    receptor_1 <- receptor_2 <- secreted_1 <- secreted_2 <- L_3 <- LR_ID <- NULL
  data <- LRcp_raw
  dt_interac <- setDT(data$interaction_table)
  dt_multi <- setDT(data$multidata_table)
  dt_prot <- setDT(data$protein_table)
  dt_gene <- setDT(data$gene_table)
  dt_comp <- setDT(data$complex_composition_table)
  dt_full <- merge.data.table(
    x = dt_interac,
    y = dt_multi,
    by.x = "multidata_1_id",
    by.y = "id_multidata",
    all.x = TRUE,
    sort = FALSE
  )
  dt_full <- merge.data.table(
    x = dt_full,
    y = dt_multi,
    by.x = "multidata_2_id",
    by.y = "id_multidata",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_1", "_2")
  )
  dt_full <- merge.data.table(
    x = dt_full,
    y = dt_prot,
    by.x = "multidata_1_id",
    by.y = "protein_multidata_id",
    all.x = TRUE,
    sort = FALSE
  )
  dt_full <- merge.data.table(
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
  dt_full <- merge.data.table(
    x = dt_full,
    y = dt_comp_dc,
    by.x = "multidata_1_id",
    by.y = "complex_multidata_id",
    all.x = TRUE,
    sort = FALSE
  )
  dt_full <- merge.data.table(
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
  dt_gene$temp_id <- rowid(dt_gene$protein_id)
  dt_gene <- dt_gene[temp_id == 1, ]
  dt_gene[, temp_id := NULL]
  dt_full[, c(paste0("id_gene_multi_", 1:3, "_1"), paste0("id_gene_multi_", 1:3, "_2")) := c(
    lapply(1:3, function(i) {
      dt_gene[.SD, on = paste0("protein_id==id_protein_multi_", i, "_1"), get("x.hgnc_symbol")]
    }),
    lapply(1:3, function(i) {
      dt_gene[.SD, on = paste0("protein_id==id_protein_multi_", i, "_2"), get("x.hgnc_symbol")]
    })
  )]
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
  if (deconvoluted) {
    dt_full <- do.call("rbind", lapply(1:2, function(i) {
      do.call("rbind", lapply(1:3, function(j) {
        cols <- c(paste0("L_", i), paste0("R_", j))
        df <- stats::na.omit(dt_full[, cols, with = FALSE])
        colnames(df) <- c("L_1", "R_1")
        df$LR_ID <- paste(df$L_1, df$R_1, sep = "_")
        df <- df[!duplicated(df$LR_ID), ]
      }))
    }))
    dt_full <- dt_full[!duplicated(dt_full$LR_ID), ]
    dt_full[, LR_ID := NULL]
  }
  return(dt_full)
}

prepare_LR_cpdb <- function(
  species,
  one2one = FALSE,
  deconvoluted = FALSE,
  keep_id = FALSE
) {
  LR_SORTED <- NULL
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  LR <- create_LR_cpdb(
    deconvoluted = deconvoluted
  )
  if (species == "mouse") {
    genes_temp <- unique(unlist(LR))
    genes_temp <- genes_temp[!is.na(genes_temp)]
    ortho <- get_orthologs(
      genes = genes_temp,
      input_species = "human",
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    if (deconvoluted) {
      nL <- 1
      nR <- 1
      LR$SOURCE <- "CPDB_DECONV"
      cols_to_keep <- c(
        "LR_SORTED", "SOURCE",
        "LIGAND_1", "RECEPTOR_1",
        "LIGAND_1_CONF", "RECEPTOR_1_CONF",
        "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
      )
    } else {
      nL <- 2
      nR <- 3
      LR$SOURCE <- "CPDB"
      cols_to_keep <- c(
        "LR_SORTED", "SOURCE",
        "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
        "LIGAND_1_CONF", "LIGAND_2_CONF", "RECEPTOR_1_CONF", "RECEPTOR_2_CONF", "RECEPTOR_3_CONF",
        "LIGAND_1_TYPE", "LIGAND_2_TYPE", "RECEPTOR_1_TYPE", "RECEPTOR_2_TYPE", "RECEPTOR_3_TYPE"
      )
      if (keep_id) {
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
  }
  if (species == "human") {
    setnames(
      x = LR,
      old = c(
        "L_1", "L_2",
        "R_1", "R_2", "R_3"
      ),
      new = c(
        "LIGAND_1", "LIGAND_2",
        "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
      )
    )
    LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
      temp <- c(
        sapply(1:2, function(j) {
          get(paste0("LIGAND_", j))[[i]]
        }),
        sapply(1:3, function(j) {
          get(paste0("RECEPTOR_", j))[[i]]
        })
      )
      temp <- temp[!is.na(temp)]
      temp <- sort(temp)
      temp <- paste0(temp, collapse = "_")
    }))]
    LR <- LR[!duplicated(LR_SORTED)]
    LR$SOURCE <- "CPDB"
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "LIGAND_2",
      "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
    )
  }
  return(LR[, cols_to_keep, with = FALSE])
}

prepare_LR_CellChat <- function(
  species
) {
  SOURCE <- NULL
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("Package \"CellChat\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  LIGAND_1 <- RECEPTOR_1 <- RECEPTOR_2 <- LR_SORTED <- interaction_name_2 <- temp <- new <- NULL
  if (species == "mouse") {
    LR <- CellChat::CellChatDB.mouse$interaction
  }
  if (species == "human") {
    LR <- CellChat::CellChatDB.human$interaction
  }
  setDT(LR)
  setnames(
    x = LR,
    old = c("evidence", "annotation"),
    new = c("SOURCE", "ANNOTATION")
  )
  LR[, LIGAND_1 := sub(" - .*", "", interaction_name_2)]
  LR[, temp := sub(".* - ", "", interaction_name_2)]
  LR[, RECEPTOR_1 := ifelse(grepl("+", temp, fixed = TRUE), gsub(".*\\((.+)\\+.*", "\\1", temp), temp)]
  LR[, RECEPTOR_2 := ifelse(grepl("+", temp, fixed = TRUE), gsub(".*\\+(.+)\\).*", "\\1", temp), NA)]
  LR[, temp := NULL]
  LR[, LIGAND_1 := gsub(" ", "", LIGAND_1)]
  LR[, RECEPTOR_1 := gsub(" ", "", RECEPTOR_1)]
  LR[, RECEPTOR_2 := gsub(" ", "", RECEPTOR_2)]
  if (species == "mouse") {
    # some CELLCHAT gene names (70) are not mgi_symbols and we need to convert them manually...
    convert_table <- CellChat_conversion
    genes_to_rm <- convert_table[new == "remove"]
    genes_to_change <- convert_table[new != "remove"]
    LR <- LR[!(LIGAND_1 %in% genes_to_rm$old) & !(RECEPTOR_1 %in% genes_to_rm$old) & !(RECEPTOR_2 %in% genes_to_rm$old)]
    LR[genes_to_change,
       `:=`(LIGAND_1 = new),
       on = "LIGAND_1==old"
       ][
         genes_to_change,
         `:=`(RECEPTOR_1 = new),
         on = "RECEPTOR_1==old"
         ][
           genes_to_change,
           `:=`(RECEPTOR_2 = new),
           on = "RECEPTOR_2==old"
           ]
  }
  LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
    temp <- c(LIGAND_1[[i]], RECEPTOR_1[[i]], RECEPTOR_2[[i]])
    temp <- temp[!is.na(temp)]
    temp <- sort(temp)
    temp <- paste0(temp, collapse = "_")
  }))]
  LR <- LR[!duplicated(LR_SORTED)]
  cols_to_keep <- c(
    "LR_SORTED",
    "ANNOTATION", "SOURCE",
    "LIGAND_1", "RECEPTOR_1", "RECEPTOR_2"
  )
  LR <- LR[, cols_to_keep, with = FALSE]
  LR[, SOURCE := gsub(" ", "", SOURCE)]
  return(LR)
}

prepare_LR_ICELLNET <- function(
  species,
  one2one = FALSE
) {
  LR_SORTED <- SOURCE <- R3 <- NULL
  if (!(species %in% c("human", "mouse"))) {
    stop("`species` muste be either 'mouse' or 'human'")
  }
  if (!requireNamespace("curl", quietly = TRUE)) {
    stop("Package \"curl\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  LR <- utils::read.csv(
    file = curl::curl(url = "https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/database.tsv"),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = ""
  )
  setDT(LR)
  setnames(
    x = LR,
    old = c(
      "Ligand 1", "Ligand 2",
      "Receptor 1", "Receptor 2", "Receptor 3",
      "PubMed ID", "Family", "Subfamily", "Classifications"
    ),
    new = c(
      "L1", "L2",
      "R1", "R2", "R3",
      "SOURCE", "FAMILY", "SUBFAMILY", "ANNOTATION"
    )
  )
  LR[, R3 := ifelse(R3 %in% c(" ", "    "), NA, R3)]
  if (species == "mouse") {
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
  }
  if (species == "human") {
    setnames(
      x = LR,
      old = c(
        "L1", "L2",
        "R1", "R2", "R3"
      ),
      new = c(
        "LIGAND_1", "LIGAND_2",
        "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
      )
    )
    LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
      temp <- c(
        sapply(1:2, function(j) {
          get(paste0("LIGAND_", j))[[i]]
        }),
        sapply(1:3, function(j) {
          get(paste0("RECEPTOR_", j))[[i]]
        })
      )
      temp <- temp[!is.na(temp)]
      temp <- sort(temp)
      temp <- paste0(temp, collapse = "_")
    }))]
    LR <- LR[!duplicated(LR_SORTED)]
    cols_to_keep <- c(
      "LR_SORTED",
      "FAMILY", "SUBFAMILY", "ANNOTATION", "SOURCE",
      "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
    )
  }
  LR <- LR[!is.na(SOURCE)]
  LR[, SOURCE := gsub(" ", "", SOURCE)]
  LR[, SOURCE := paste0("PMID:", SOURCE)]
  LR[, SOURCE := gsub(";", ";PMID:", SOURCE)]
  return(LR[, cols_to_keep, with = FALSE])
}

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
    sapply(1:nL, function(i) {
      paste0("LIGAND_", i, c("", "_CONF", "_TYPE"))
    }),
    sapply(1:nR, function(i) {
      paste0("RECEPTOR_", i, c("", "_CONF", "_TYPE"))
    })
  )
  merge_id <- c("mouse_symbol", "confidence", "type")
  LR_temp <- copy(LR_dt)
  LR_temp[, c(out_names) :=
            c(
              sapply(
                1:nL,
                function(i) {
                  as.list(
                    ortho_dt[.SD,
                             on = c(paste0("human_symbol==", charL, i)),
                             mget(paste0("x.", merge_id))
                             ]
                  )
                }
              ),
              sapply(
                1:nR,
                function(i) {
                  as.list(
                    ortho_dt[.SD,
                             on = c(paste0("human_symbol==", charR, i)),
                             mget(paste0("x.", merge_id))
                             ]
                  )
                }
              )
            )]
  LR_temp <- stats::na.omit(LR_temp, cols = c("LIGAND_1", "RECEPTOR_1"))
  LR_temp[, to_keep := sapply(1:nrow(.SD), function(i) {
    all(c(
      sapply(1:nL, function(j) {
        !(is.na(get(paste0("LIGAND_", j))[[i]]) & !is.na(get(paste0(charL, j))[[i]]))
      }),
      sapply(1:nR, function(j) {
        !(is.na(get(paste0("RECEPTOR_", j))[[i]]) & !is.na(get(paste0(charR, j))[[i]]))
      })
    ))
  })]
  LR_temp <- LR_temp[to_keep == TRUE]
  LR_temp[, to_keep := NULL]
  LR_temp[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
    temp <- c(
      sapply(1:nL, function(j) {
        get(paste0("LIGAND_", j))[[i]]
      }),
      sapply(1:nR, function(j) {
        get(paste0("RECEPTOR_", j))[[i]]
      })
    )
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
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package \"biomaRt\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  ensembl_gene_id <- inl <- outl <- output <- input <- confidence <- NULL
  if (input_species == "mouse") {
    id_in <- "mmusculus"
    id_out <- "hsapiens"
    id_gene <- "mgi_symbol"
    name_in <- "mouse"
    name_out <- "human"
  } else if (input_species == "human") {
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
    value = ensembl$ensembl_gene_id
  )
  setDT(ensembl)
  setDT(ensembl_conv)
  ensembl_all <- merge.data.table(
    x = ensembl,
    y = ensembl_conv,
    by = "ensembl_gene_id",
    all = TRUE,
    sort = FALSE
  )
  if (one2one) {
    ensembl_all <- ensembl_all[eval(as.symbol(ortho_type)) == "ortholog_one2one", ]
  } else {
    # ensembl_all <- ensembl_all[eval(as.symbol(ortho_type)) %in% c('ortholog_one2one','ortholog_one2many'),]
  }
  ensembl_all <- ensembl_all[, ensembl_gene_id := NULL]
  ensembl_all <- unique(ensembl_all)
  setnames(
    x = ensembl_all,
    old = c(id_gene, gene_name, ortho_confidence, ortho_type),
    new = c("input", "output", "confidence", "type")
  )
  ensembl_all <- stats::na.omit(ensembl_all)
  # check for remaining duplicate
  if (sum(duplicated(ensembl_all[["input"]])) > 0) {
    ensembl_all[, inl := tolower(input)]
    ensembl_all[, outl := tolower(output)]
    dup_input <- unique(ensembl_all$input[duplicated(ensembl_all$input)])
    for (g in dup_input) {
      dt_g <- ensembl_all[input == g]
      dt_gconf <- dt_g[confidence == 1]
      if (nrow(dt_gconf) == 1) {
        ensembl_all <- ensembl_all[!(input == g & confidence == 0)]
      } else if (nrow(dt_gconf) == 0) {
        dt_gsame <- dt_g[inl == outl]
        if (nrow(dt_gsame) == 0) {
          g_keep <- dt_g[1]$output
        } else {
          g_keep <- dt_gsame[1]$output
        }
        ensembl_all <- ensembl_all[!(input == g & output != g_keep)]
      } else {
        dt_gsame <- dt_gconf[inl == outl]
        if (nrow(dt_gsame) == 0) {
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
    if (sum(duplicated(ensembl_all[["input"]])) > 0) {
      warning("There are some duplicates from orthology conversion. Removing them by using 'unique'.")
      ensembl_all <- unique(ensembl_all, by = "input")
    }
  }
  setnames(
    x = ensembl_all,
    old = c("input", "output"),
    new = c(paste0(name_in, "_symbol"), paste0(name_out, "_symbol"))
  )
  return(ensembl_all)
}
