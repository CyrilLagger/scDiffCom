build_LRI <- function(
  species = c("mouse", "human", "rat")
) {
  species <- match.arg(species)
  LRI_not_curated <- combine_LR_db(
    species = species,
    one2one = FALSE,
    curated = FALSE
  )
  LRI_curated <- combine_LR_db(
    species = species,
    one2one = FALSE,
    curated = TRUE
  )
  LRI_GO <- get_GO_interactions(
    species = species,
    LR_db = LRI_curated$LR_full
  )
  LRI_KEGG <- get_KEGG_PWS_interactions(
    species = species,
    LR_db = LRI_curated$LR_full
  )
  return(
    list(
      #LRI_not_curated = LRI_not_curated$LR_full,
      LRI_curated = LRI_curated$LR_full,
      LRI_curated_GO = LRI_GO[, c("LRI", "GO_ID")],
      LRI_curated_KEGG = LRI_KEGG,
      LRI_retrieved_dates = LRI_curated$LR_retrieved_dates,
      LRI_retrieved_from = LRI_curated$LR_retrieved_from,
      LRI_biomart_ensembl_version = "https://nov2020.archive.ensembl.org"
    )
  )
}

combine_LR_db <- function(
  species,
  one2one = FALSE,
  curated = TRUE
) {
  DATABASE <- SOURCE <- ANNOTATION <- FAMILY <- SUBFAMILY <- keep_subLR <-
    SOURCE_CLEAN <- SOURCE_no_digit <- is_complex_temp <- LIGAND_2 <-
    RECEPTOR_2 <- LR_vectorized_temp <- N_IS_SUBPART <- i.N_IS_SUBPART <-
    LRI <- LIGAND_1 <- RECEPTOR_1 <- RECEPTOR_3 <- LIGAND_1_CONF <-
    LIGAND_2_CONF <- RECEPTOR_1_CONF <- RECEPTOR_2_CONF <- RECEPTOR_3_CONF <- NULL
  # fully curated
  LR_connectomeDB2020 <- prepare_LR_connectomeDB2020(
    species = species,
    one2one = one2one
  )
  LR_CellTalkDB <- prepare_LR_CellTalkDB(
    species = species,
    one2one = one2one
  )
  LR_ICELLNET <- prepare_LR_ICELLNET(
    species = species,
    one2one = one2one
  )
  LR_CellChat <- prepare_LR_CellChat(
    species = species,
    one2one = one2one
  )
  LR_CellPhoneDB <- prepare_LR_CellPhoneDB(
    species = species,
    one2one = one2one,
    deconvoluted = FALSE,
    keep_id = FALSE
  )
  # not fully curated
  LR_SingleCellSignalR <- prepare_LR_SingleCellSignalR(
    species = species,
    one2one = one2one
  )
  LR_NicheNet <- prepare_LR_NicheNet(
    species = species,
    one2one = one2one
  )
  # LR_scTensor <- prepare_LR_scTensor(
  #   species = species
  # )
  if (curated) {
    LR_SingleCellSignalR$LR <- LR_SingleCellSignalR$LR[SOURCE != "PPI"]
    LR_NicheNet$LR <- LR_NicheNet$LR[SOURCE != "PPI"]
    # LR_scTensor$LR <- LR_scTensor$LR[SOURCE != "PPI"]
  }
  LR_retrieved_dates <- list(
    "CellChat" = LR_CellChat$retrieved_date,
    "CellPhoneDB" = LR_CellPhoneDB$retrieved_date,
    "CellTalkDB" = LR_CellTalkDB$retrieved_date,
    "connectomeDB2020" = LR_connectomeDB2020$retrieved_date,
    "ICELLNET" = LR_ICELLNET$retrieved_date,
    "NicheNet" = LR_NicheNet$retrieved_date,
    "SingleCellSignalR" = LR_SingleCellSignalR$retrieved_date#,
    #"scTensor" = LR_scTensor$retrieved_date
  )
  LR_retrieved_from <- list(
    "CellChat" = LR_CellChat$retrieved_from,
    "CellPhoneDB" = LR_CellPhoneDB$retrieved_from,
    "CellTalkDB" = LR_CellTalkDB$retrieved_from,
    "connectomeDB2020" = LR_connectomeDB2020$retrieved_from,
    "ICELLNET" = LR_ICELLNET$retrieved_from,
    "NicheNet" = LR_NicheNet$retrieved_from,
    "SingleCellSignalR" = LR_SingleCellSignalR$retrieved_from#,
    #"scTensor" = LR_scTensor$retrieved_from
  )
  LR_full <- rbindlist(
    list(
      "connectomeDB2020" = LR_connectomeDB2020$LR,
      "CellTalkDB" = LR_CellTalkDB$LR,
      #"scTensor" = LR_scTensor$LR,
      "SingleCellSignalR" = LR_SingleCellSignalR$LR,
      "NicheNet" = LR_NicheNet$LR,
      "CellPhoneDB" = LR_CellPhoneDB$LR,
      "CellChat" = LR_CellChat$LR,
      "ICELLNET" = LR_ICELLNET$LR
    ),
    use.names = TRUE,
    fill = TRUE,
    idcol = "DATABASE"
  )
  col_md <- colnames(LR_full)
  col_md <- col_md[col_md != "LR_SORTED"]
  db_names <- c(
    "connectomeDB2020",
    "CellTalkDB",
    "CellChat",
    "CellPhoneDB",
    "SingleCellSignalR",
    #"scTensor",
    "ICELLNET",
    "NicheNet"
  )
  LR_full <- dcast.data.table(
    LR_full,
    formula = LR_SORTED ~ DATABASE,
    value.var = col_md
  )
  LR_full[
    ,
    paste0("LIGAND_", 1:2) := lapply(
      1:2,
      function(i) {
        ifelse(
          !is.na(get(paste0("LIGAND_", i, "_CellPhoneDB"))),
          get(paste0("LIGAND_", i, "_CellPhoneDB")),
          ifelse(
            !is.na(get(paste0("LIGAND_", i, "_ICELLNET"))),
            get(paste0("LIGAND_", i, "_ICELLNET")),
            ifelse(
              !is.na(get(paste0("LIGAND_", i, "_SingleCellSignalR"))),
              get(paste0("LIGAND_", i, "_SingleCellSignalR")),
              ifelse(
                #!is.na(get(paste0("LIGAND_", i, "_scTensor"))),
                #get(paste0("LIGAND_", i, "_scTensor")),
                #ifelse(
                !is.na(get(paste0("LIGAND_", i, "_NicheNet"))),
                get(paste0("LIGAND_", i, "_NicheNet")),
                ifelse(
                  !is.na(get(paste0("LIGAND_", i, "_connectomeDB2020"))),
                  get(paste0("LIGAND_", i, "_connectomeDB2020")),
                  ifelse(
                    !is.na(get(paste0("LIGAND_", i, "_CellTalkDB"))),
                    get(paste0("LIGAND_", i, "_CellTalkDB")),
                    get(paste0("LIGAND_", i, "_CellChat"))
                  )
                )
                #)
              )
            )
          )
        )
      }
    )
  ]
  LR_full[
    ,
    paste0("RECEPTOR_", 1:3) := lapply(
      1:3,
      function(i) {
        ifelse(
          !is.na(get(paste0("RECEPTOR_", i, "_CellPhoneDB"))),
          get(paste0("RECEPTOR_", i, "_CellPhoneDB")),
          ifelse(
            !is.na(get(paste0("RECEPTOR_", i, "_ICELLNET"))),
            get(paste0("RECEPTOR_", i, "_ICELLNET")),
            ifelse(
              !is.na(get(paste0("RECEPTOR_", i, "_SingleCellSignalR"))),
              get(paste0("RECEPTOR_", i, "_SingleCellSignalR")),
              ifelse(
                #!is.na(get(paste0("RECEPTOR_", i, "_scTensor"))),
                #get(paste0("RECEPTOR_", i, "_scTensor")),
                #ifelse(
                !is.na(get(paste0("RECEPTOR_", i, "_NicheNet"))),
                get(paste0("RECEPTOR_", i, "_NicheNet")),
                ifelse(
                  !is.na(get(paste0("RECEPTOR_", i, "_connectomeDB2020"))),
                  get(paste0("RECEPTOR_", i, "_connectomeDB2020")),
                  ifelse(
                    !is.na(get(paste0("RECEPTOR_", i, "_CellTalkDB"))),
                    get(paste0("RECEPTOR_", i, "_CellTalkDB")),
                    get(paste0("RECEPTOR_", i, "_CellChat"))
                  )
                )
                #)
              )
            )
          )
        )
      }
    )
  ]
  col_database <- paste0("DATABASE.1_", db_names)
  LR_full[
    ,
    c("DATABASE") := do.call(paste, c(.SD, sep = ";")),
    .SDcols = col_database
  ]
  LR_full[, c("DATABASE") := gsub("NA|NA;|;NA", "", DATABASE)]
  LR_full[, c("N_DB") := rowSums(!is.na(.SD)), .SDcols = col_database]
  col_source <- paste0("SOURCE_", db_names)
  LR_full[
    ,
    c("SOURCE") := do.call(paste, c(.SD, sep = ";")),
    .SDcols = col_source
  ]
  LR_full[, c("SOURCE") := gsub("NA|NA;|;NA", "", SOURCE)]
  col_anno <- paste0("ANNOTATION_", db_names)
  LR_full[
    ,
    c("ANNOTATION") := do.call(paste, c(.SD, sep = ";")),
    .SDcols = col_anno
  ]
  LR_full[, c("ANNOTATION") := gsub("NA|NA;|;NA", "", ANNOTATION)]
  col_fam <- paste0("FAMILY_", db_names)
  LR_full[, c("FAMILY") := do.call(paste, c(.SD, sep = ";")), .SDcols = col_fam]
  LR_full[, c("FAMILY") := gsub("NA|NA;|;NA", "", FAMILY)]
  col_subfam <- paste0("SUBFAMILY_", db_names)
  LR_full[
    , c("SUBFAMILY") := do.call(paste, c(.SD, sep = ";")),
    .SDcols = col_subfam]
  LR_full[, c("SUBFAMILY") := gsub("NA|NA;|;NA", "", SUBFAMILY)]
  if(species %in% c("mouse", "rat")) {
    for (id_loop in c(paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3))) {
      cols_conf <- paste0(id_loop, "_CONF_", db_names)
      LR_full[
        ,
        paste0(id_loop, "_CONF") := do.call(pmax, c(.SD, na.rm = TRUE)),
        .SDcols = cols_conf]
      LR_full[
        , paste0(id_loop, "_CONF") := ifelse(
          is.na(get(paste0(id_loop, "_CONF"))) & !is.na(get(id_loop)),
          1,
          get(paste0(id_loop, "_CONF"))
        )]
      cols_type <- paste0(id_loop, "_TYPE_", db_names)
      LR_full[, paste0(id_loop, "_TYPE") := do.call(
        pmax,
        c(.SD, na.rm = TRUE)
      ),
      .SDcols = cols_type]
      LR_full[
        ,
        paste0(id_loop, "_TYPE") := ifelse(
          is.na(get(paste0(id_loop, "_TYPE"))) & !is.na(get(id_loop)),
          "provided",
          get(paste0(id_loop, "_TYPE")
          )
        )
      ]
    }
  }
  if (curated) {
    LR_full[
      ,
      is_complex_temp := fifelse(
        !is.na(LIGAND_2) | !is.na(RECEPTOR_2),
        TRUE,
        FALSE
      )
    ]
    LR_full[
      ,
      LR_vectorized_temp := list(
        sapply(
          1:nrow(.SD),
          function(i) {
            temp <- c(
              sapply(
                1:2,
                function(j) {
                  get(paste0("LIGAND_", j))[[i]]
                }
              ),
              sapply(
                1:3, function(j) {
                  get(paste0("RECEPTOR_", j))[[i]]
                }
              )
            )
            temp <- temp[!is.na(temp)]
          }
        )
      )
    ]
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
    LR_full[, keep_subLR := grepl("CellPhoneDB|CellChat|ICELLNET", DATABASE)]
    LR_full <- LR_full[N_IS_SUBPART == 0 | keep_subLR == TRUE]
    #LR_full <- LR_full[DATABASE != "scTensor"]
    cols_to_keep <- NULL
  } else {
    cols_to_keep <- NULL
  }
  LR_full[
    ,
    LRI := list(
      sapply(
        1:nrow(.SD),
        function(i) {
          temp1 <- c(LIGAND_1[[i]], LIGAND_2[[i]])
          temp1 <- temp1[!is.na(temp1)]
          temp1 <- paste0(temp1, collapse = "_")
          temp2 <- c(RECEPTOR_1[[i]], RECEPTOR_2[[i]], RECEPTOR_3[[i]])
          temp2 <- temp2[!is.na(temp2)]
          temp2 <- paste0(temp2, collapse = "_")
          return(paste(temp1, temp2, sep = ":"))
        }
      )
    )
  ]
  if (species %in% c("mouse", "rat")) {
    if (curated) {
      # remove LRI with 0 orthology confidence
      # if not provided by another database
      genes_conf <- unique(
        c(
          LR_full[LIGAND_1_CONF == 1]$LIGAND_1,
          LR_full[LIGAND_2_CONF == 1]$LIGAND_2,
          LR_full[RECEPTOR_1_CONF == 1]$RECEPTOR_1,
          LR_full[RECEPTOR_2_CONF == 1]$RECEPTOR_2,
          LR_full[RECEPTOR_3_CONF == 1]$RECEPTOR_3
        )
      )
      cols_conf <- c(
        "LIGAND_1", "LIGAND_2",
        "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
      )
      LR_full[
        ,
        paste0(cols_conf, "_CONF") := lapply(
          cols_conf,
          function(i) {
            ifelse(
              get(i) %in% genes_conf, 1, get(paste0(i, "_CONF"))
            )
          }
        )
      ]
      LR_full <- LR_full[
        LIGAND_1_CONF == 1 &
          (LIGAND_2_CONF == 1 | is.na(LIGAND_2_CONF)) &
          RECEPTOR_1_CONF == 1 &
          (RECEPTOR_2_CONF == 1 | is.na (RECEPTOR_2_CONF)) &
          (RECEPTOR_3_CONF == 1 | is.na (RECEPTOR_3_CONF))
      ]
    }
    cols_to_keep <- c(
      cols_to_keep,
      "LRI",
      paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3),
      "DATABASE", "SOURCE"#,
      #paste0("LIGAND_", 1:2, "_CONF"), paste0("LIGAND_", 1:2, "_TYPE"),
      #paste0("RECEPTOR_", 1:3, "_CONF"), paste0("RECEPTOR_", 1:3, "_TYPE")
    )
  }
  if (species == "human") {
    cols_to_keep <- c(
      cols_to_keep,
      "LRI",
      paste0("LIGAND_", 1:2), paste0("RECEPTOR_", 1:3),
      "DATABASE", "SOURCE"
    )
  }
  LR_full <- LR_full[, cols_to_keep, with = FALSE]
  #clean the columns DATABASE and SOURCE
  LR_full[
    ,
    DATABASE := sapply(
      1:nrow(.SD),
      function(i) {
        temp <- unlist(strsplit(.SD[i,]$DATABASE, ";"))
        paste0(sort(unique(temp)), collapse = ";")
      }
    )
  ]
  LR_full[, SOURCE := gsub("\u00a0", "", SOURCE)]
  LR_full[, SOURCE := gsub(")", "", SOURCE, fixed = TRUE)]
  LR_full[, SOURCE := gsub(";;", ";", SOURCE, fixed = TRUE)]
  LR_full[
    ,
    SOURCE := sapply(
      1:nrow(.SD),
      function(i) {
        temp <- unlist(strsplit(.SD[i,]$SOURCE, ";"))
        paste0(sort(unique(temp)), collapse = ";")
      }
    )
  ]
  return(
    list(
      LR_full = LR_full,
      LR_retrieved_dates = LR_retrieved_dates,
      LR_retrieved_from = LR_retrieved_from
    )
  )
}

prepare_LR_connectomeDB2020 <- function(
  species,
  one2one
) {
  LR_SORTED <- SOURCE <- L1 <- R1 <- NULL
  #Last time updated
  retrieved_date <- as.Date("2021-03-22")
  retrieved_from <- "https://asrhou.github.io/NATMI/"
  #we checked manually that they are all HGCN Approved Symbol
  LR <- NATMI_connectomeDB2020
  setDT(LR)
  setnames(
    x = LR,
    old = c("Ligand.gene.symbol", "Receptor.gene.symbol", "PMID.support"),
    new = c("L1", "R1", "SOURCE")
  )
  #remove genes that should be excluded according to our manual curation
  LR <- LR[
    !(L1 %in% GENES_to_remove_human$gene) &
      !(R1 %in% GENES_to_remove_human$gene)
  ]
  if (species %in% c("mouse", "rat")) {
    ortho <- get_orthologs(
      genes = unique(c(LR$L1, LR$R1)),
      input_species = "human",
      output_species = species,
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    LR <- merge_LR_orthologs(
      LR_dt = LR,
      ortho_dt = ortho,
      nL = 1,
      charL = "L",
      nR = 1,
      charR = "R",
      output_species = species,
      input_species = "human"
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
  return(
    list(
      retrieved_date = retrieved_date,
      retrieved_from = retrieved_from,
      LR = LR[, cols_to_keep, with = FALSE]
    )
  )
}

prepare_LR_CellTalkDB <- function(
  species,
  one2one
) {
  LR_SORTED <- SOURCE <- LIGAND_1 <- RECEPTOR_1 <- NULL
  retrieved_date <- as.Date("2021-03-22")
  retrieved_from <- "http://tcm.zju.edu.cn/celltalkdb/"
  if (species %in% c("mouse", "rat")) {
    LR <- CellTalkDB_mouse
  }
  if (species == "human") {
    LR <- CellTalkDB_human
    #we checked manually that they are all HGCN Approved Symbol
  }
  setDT(LR)
  setnames(
    LR,
    old = c("ligand_gene_symbol", "receptor_gene_symbol", "evidence"),
    new = c("LIGAND_1", "RECEPTOR_1", "SOURCE")
  )
  #remove genes that should be excluded according to our manual curation
  if (species %in% c("mouse", "rat")) {
    LR <- LR[
      !(LIGAND_1 %in% GENES_to_remove_mouse$gene) &
        !(RECEPTOR_1 %in% GENES_to_remove_mouse$gene)
    ]
    if (species == "mouse") {
      #there are genes that are not MGI approved symbol, we change them
      LR[LR == "Il1f8"] <- "Il36b"
      LR[LR == "Il1f6"] <- "Il36a"
      LR[LR == "Il1f5"] <- "Il36rn"
    }
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1"
    )
  }
  if (species == "human") {
    LR <- LR[
      !(LIGAND_1 %in% GENES_to_remove_human$gene) &
        !(RECEPTOR_1 %in% GENES_to_remove_human$gene)
    ]
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1"
    )
  }
  if (species == "rat") {
    ortho <- get_orthologs(
      genes = unique(c(LR$LIGAND_1, LR$RECEPTOR_1)),
      input_species = "mouse",
      output_species = species,
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    LR <- merge_LR_orthologs(
      LR_dt = LR,
      ortho_dt = ortho,
      nL = 1,
      charL = "LIGAND_",
      nR = 1,
      charR = "RECEPTOR_",
      output_species = species,
      input_species = "mouse"
    )
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "RECEPTOR_1",
      "LIGAND_1_CONF", "RECEPTOR_1_CONF",
      "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
    )
  }
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
  #some manual fixes to correct bad formatting in their original DB
  LR[, SOURCE := sub(",$", "", SOURCE)]
  LR[, SOURCE := gsub(";NO", "", SOURCE, fixed = TRUE)]
  LR[, SOURCE := gsub(";", ",", SOURCE, fixed = TRUE)]
  LR[, SOURCE := gsub(",,", ",", SOURCE, fixed = TRUE)]
  LR[, SOURCE := paste0("PMID:", SOURCE)]
  LR[, SOURCE := gsub(",", ";PMID:", SOURCE, fixed = TRUE)]
  return(
    list(
      retrieved_date = retrieved_date,
      retrieved_from = retrieved_from,
      LR = LR[, cols_to_keep, with = FALSE]
    )
  )
}

prepare_LR_ICELLNET <- function(
  species,
  one2one
) {
  LR_SORTED <- SOURCE <- L1 <- L2 <- R1 <- R2 <- R3 <- NULL
  retrieved_date <- as.Date("2021-03-22")
  retrieved_from <- paste0(
    "https://raw.githubusercontent.com/soumelis-lab/",
    "ICELLNET/master/data/ICELLNETdb.tsv"
  )
  # LR <- utils::read.csv(
  #   file = curl::curl(url =  paste0(
  #     "https://raw.githubusercontent.com/soumelis-lab/",
  #     "ICELLNET/master/data/ICELLNETdb.tsv"
  #   )),
  #   sep = "\t",
  #   header = TRUE,
  #   check.names = FALSE,
  #   stringsAsFactors = FALSE,
  #   na.strings = ""
  # )
  # actually saved as internal data
  # there were two symbols that were not HGNC, probably typos in ICELLNET:
  # 3,00 NPR and VCTN1. we fixed it in ICELLNET_human
  LR <- ICELLNET_human
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
  #remove genes that should be excluded according to our manual curation
  LR <- LR[
    !(L1 %in% GENES_to_remove_human$gene) &
      !(L2 %in% GENES_to_remove_human$gene) &
      !(R1 %in% GENES_to_remove_human$gene) &
      !(R2 %in% GENES_to_remove_human$gene) &
      !(R3 %in% GENES_to_remove_human$gene)
  ]
  if (species %in% c("mouse", "rat")) {
    ortho <- get_orthologs(
      genes = unique(c(LR$L1, LR$L2, LR$R1, LR$R2, LR$R3)),
      input_species = "human",
      output_species = species,
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    LR <- merge_LR_orthologs(
      LR_dt = LR,
      ortho_dt = ortho,
      nL = 2,
      charL = "L",
      nR = 3,
      charR = "R",
      output_species = species,
      input_species = "human"
    )
    cols_to_keep <- c(
      "LR_SORTED",
      "FAMILY", "SUBFAMILY", "ANNOTATION", "SOURCE",
      "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
      "LIGAND_1_CONF", "LIGAND_2_CONF",
      "RECEPTOR_1_CONF", "RECEPTOR_2_CONF", "RECEPTOR_3_CONF",
      "LIGAND_1_TYPE", "LIGAND_2_TYPE",
      "RECEPTOR_1_TYPE", "RECEPTOR_2_TYPE", "RECEPTOR_3_TYPE"
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
  return(
    list(
      retrieved_date = retrieved_date,
      retrieved_from = retrieved_from,
      LR = LR[, cols_to_keep, with = FALSE]
    )
  )
}

prepare_LR_CellChat <- function(
  species,
  one2one
) {
  LIGAND_1 <- RECEPTOR_1 <- RECEPTOR_2 <- LR_SORTED <- SOURCE <-
    interaction_name_2 <- temp <- new <- NULL
  # if (!requireNamespace("CellChat", quietly = TRUE)) {
  #   stop(
  #     paste0(
  #       "Package \"CellChat\" needed for this function to work.",
  #       "Please install it."
  #     ),
  #     call. = FALSE
  #   )
  # }
  retrieved_date <- as.Date("2021-03-22")
  retrieved_from <- "CellChat::CellChatDB"
  if (species %in% c("mouse", "rat")) {
    #LR <- CellChat::CellChatDB.mouse$interaction
    LR <- CellChat_mouse
  }
  if (species == "human") {
    #LR <- CellChat::CellChatDB.human$interaction
    LR <- CellChat_human
  }
  setDT(LR)
  setnames(
    x = LR,
    old = c("evidence", "annotation"),
    new = c("SOURCE", "ANNOTATION")
  )
  LR[, LIGAND_1 := sub(" - .*", "", interaction_name_2)]
  LR[, temp := sub(".* - ", "", interaction_name_2)]
  LR[, RECEPTOR_1 := ifelse(grepl("+", temp, fixed = TRUE),
                            gsub(".*\\((.+)\\+.*", "\\1", temp), temp)]
  LR[, RECEPTOR_2 := ifelse(grepl("+", temp, fixed = TRUE),
                            gsub(".*\\+(.+)\\).*", "\\1", temp), NA)]
  LR[, temp := NULL]
  LR[, LIGAND_1 := gsub(" ", "", LIGAND_1)]
  LR[, RECEPTOR_1 := gsub(" ", "", RECEPTOR_1)]
  LR[, RECEPTOR_2 := gsub(" ", "", RECEPTOR_2)]
  # some CellChat gene names  are not mgi/HGNC_symbols
  if (species %in% c("mouse", "rat")) {
    convert_table <- CellChat_conversion_mouse
  }
  if (species == "human") {
    convert_table <- CellChat_conversion_human
  }
  genes_to_rm <- convert_table[new == "remove"]
  genes_to_change <- convert_table[new != "remove"]
  LR <- LR[!(LIGAND_1 %in% genes_to_rm$old) &
             !(RECEPTOR_1 %in% genes_to_rm$old) &
             !(RECEPTOR_2 %in% genes_to_rm$old)]
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
  #remove genes that should be excluded according to our manual curation
  if (species %in% c("mouse", "rat")) {
    LR <- LR[
      !(LIGAND_1 %in% GENES_to_remove_mouse$gene) &
        !(RECEPTOR_1 %in% GENES_to_remove_mouse$gene) &
        !(RECEPTOR_2 %in% GENES_to_remove_mouse$gene)
    ]
  }
  if (species == "human") {
    LR <- LR[
      !(LIGAND_1 %in% GENES_to_remove_human$gene) &
        !(RECEPTOR_1 %in% GENES_to_remove_human$gene) &
        !(RECEPTOR_2 %in% GENES_to_remove_human$gene)
    ]
  }
  if (species == "rat") {
    ortho <- get_orthologs(
      genes = unique(c(LR$LIGAND_1, LR$RECEPTOR_1, LR$RECEPTOR_2)),
      input_species = "mouse",
      output_species = species,
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    LR <- merge_LR_orthologs(
      LR_dt = LR,
      ortho_dt = ortho,
      nL = 1,
      charL = "LIGAND_",
      nR = 2,
      charR = "RECEPTOR_",
      output_species = species,
      input_species = "mouse"
    )
    cols_to_keep <- c(
      "LR_SORTED", "ANNOTATION", "SOURCE",
      "LIGAND_1", "RECEPTOR_1", "RECEPTOR_2",
      "LIGAND_1_CONF", "RECEPTOR_1_CONF", "RECEPTOR_2_CONF",
      "LIGAND_1_TYPE", "RECEPTOR_1_TYPE", "RECEPTOR_2_TYPE"
    )
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
  return(
    list(
      retrieved_date = retrieved_date,
      retrieved_from = retrieved_from,
      LR = LR[, cols_to_keep, with = FALSE]
    )
  )
}

prepare_LR_CellPhoneDB <- function(
  species,
  one2one,
  deconvoluted,
  keep_id
) {
  LR_SORTED <- L_1 <- L_2 <-
    R_1 <- R_2 <- R_3 <- NULL
  retrieved_date <- as.Date("2021-03-22")
  retrieved_from <- paste0(
    "https://github.com/Teichlab/cellphonedb-data/blob/master/cellphone.db"
  )
  LR <- create_LR_CellPhoneDB(
    deconvoluted = deconvoluted
  )
  #there are 3 genes that are not HGNC approved symbol, we change them
  LR[LR == "NOV"] <- "CCN3"
  LR[LR == "WISP3"] <- "CCN6"
  LR[LR == "YARS"] <- "YARS1"
  #remove genes that should be excluded according to our manual curation
  LR <- LR[
    !(L_1 %in% GENES_to_remove_human$gene) &
      !(L_2 %in% GENES_to_remove_human$gene) &
      !(R_1 %in% GENES_to_remove_human$gene) &
      !(R_2 %in% GENES_to_remove_human$gene) &
      !(R_3 %in% GENES_to_remove_human$gene)
  ]
  if (species %in% c("mouse", "rat")) {
    genes_temp <- unique(unlist(LR))
    genes_temp <- genes_temp[!is.na(genes_temp)]
    ortho <- get_orthologs(
      genes = genes_temp,
      input_species = "human",
      output_species = species,
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    if (deconvoluted) {
      nL <- 1
      nR <- 1
      LR$SOURCE <- "CellPhoneDB_DECONV"
      cols_to_keep <- c(
        "LR_SORTED", "SOURCE",
        "LIGAND_1", "RECEPTOR_1",
        "LIGAND_1_CONF", "RECEPTOR_1_CONF",
        "LIGAND_1_TYPE", "RECEPTOR_1_TYPE"
      )
    } else {
      nL <- 2
      nR <- 3
      LR$SOURCE <- "CellPhoneDB"
      cols_to_keep <- c(
        "LR_SORTED", "SOURCE",
        "LIGAND_1", "LIGAND_2", "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
        "LIGAND_1_CONF", "LIGAND_2_CONF", "RECEPTOR_1_CONF",
        "RECEPTOR_2_CONF", "RECEPTOR_3_CONF",
        "LIGAND_1_TYPE", "LIGAND_2_TYPE", "RECEPTOR_1_TYPE",
        "RECEPTOR_2_TYPE", "RECEPTOR_3_TYPE"
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
      charR = "R_",
      input_species = "human",
      output_species = species
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
    LR$SOURCE <- "CellPhoneDB"
    cols_to_keep <- c(
      "LR_SORTED", "SOURCE",
      "LIGAND_1", "LIGAND_2",
      "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3"
    )
  }
  return(
    list(
      retrieved_date = retrieved_date,
      retrieved_from = retrieved_from,
      LR = LR[, cols_to_keep, with = FALSE]
    )
  )
}

create_LR_CellPhoneDB <- function(
  deconvoluted
) {
  id_protein_multi_1_1 <- id_protein_multi_1_2 <- id_protein_1 <-
    id_protein_2 <- temp_id <- receptor_1 <- receptor_2 <- secreted_1 <-
    secreted_2 <- L_3 <- LR_ID <- NULL
  # retrieving CPDB db file and converting it to a list of data.frame
  # done externally of scDiffCom
  # cpdb_db_remote_path <- paste0(
  #   "https://github.com/Teichlab/cellphonedb-data/blob/master/cellphone.db"
  # )
  #cpdb_db_local_path <- "cpdb_db_local_path"
  #sqlite.driver <- RSQLite::dbDriver("SQLite")
  # cpdb_internal <- RSQLite::dbConnect(
  #   sqlite.driver,
  #   dbname = cpdb_db_local_path
  #   )
  # cpdb_table_names <- dbListTables(cpdb_internal)
  # CellPhoneDB_data <- sapply(
  #   cpdb_table_names,
  #   function(i) {
  #     RSQLite::dbReadTable(cpdb_internal, i)
  #   },
  #   USE.NAMES = TRUE,
  #   simplify = FALSE
  # )
  # saving CellPhoneDB_data as part of sydata.rda
  data <- CellPhoneDB_data
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
  dt_comp$id_comp <- paste0(
    "id_protein_multi_",
    rowid(dt_comp$complex_multidata_id)
  )
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
  dt_full[, id_protein_multi_1_1 := ifelse(
    is.na(id_protein_1), id_protein_multi_1_1, id_protein_1)]
  dt_full[, id_protein_multi_1_2 := ifelse(
    is.na(id_protein_2), id_protein_multi_1_2, id_protein_2)]
  dt_gene <- stats::na.omit(unique(dt_gene[, c("protein_id", "hgnc_symbol")]))
  dt_gene$temp_id <- rowid(dt_gene$protein_id)
  dt_gene <- dt_gene[temp_id == 1, ]
  dt_gene[, temp_id := NULL]
  dt_full[
    ,
    c(
      paste0("id_gene_multi_", 1:3, "_1"),
      paste0("id_gene_multi_", 1:3, "_2")) := c(
        lapply(1:3, function(i) {
          dt_gene[
            .SD,
            on = paste0("protein_id==id_protein_multi_", i, "_1"),
            get("x.hgnc_symbol")]
        }),
        lapply(1:3, function(i) {
          dt_gene[
            .SD,
            on = paste0("protein_id==id_protein_multi_", i, "_2"),
            get("x.hgnc_symbol")]
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
  dt_full <- dt_full[
    ,
    c("id_cp_interaction", paste0("L_", 1:3), paste0("R_", 1:3))]
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

prepare_LR_SingleCellSignalR <- function(
  species,
  one2one
) {
  LR_SORTED <- SOURCE <- PMIDs <- L1 <- R1 <- NULL
  if (!requireNamespace("SingleCellSignalR", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"SingleCellSignalR\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  retrieved_date <- as.Date("2021-03-22")
  retrieved_from <- "SingleCellSignalR::LRdb"
  LR <- SingleCellSignalR::LRdb
  setDT(LR)
  #there are 19 genes that are not HGNC approved symbol, we change them
  LR[LR == "BY55"] <- "CD160"
  LR[LR == "C14orf1"] <- "ERG28"
  LR[LR == "CTGF"] <- "CCN2"
  LR[LR == "CYR61"] <- "CCN1"
  LR[LR == "HFE2"] <- "HJV"
  LR[LR == "HLAE"] <- "HLA-E"
  LR[LR == "IL8"] <- "CXCL8"
  LR[LR == "MFI2"] <- "MELTF"
  LR[LR == "MLLT4"] <- "AFDN"
  LR[LR == "MTf"] <- "MELTF"
  LR[LR == "NGFRAP1"] <- "BEX3"
  LR[LR == "NOV"] <- "CCN3"
  LR[LR == "PVRL2"] <- "NECTIN2"
  LR[LR == "RAP"] <- "LRPAP1"
  LR[LR == "SHP1"] <- "PTPN6"
  LR[LR == "SHP2"] <- "PTPN11"
  LR[LR == "TIM-1"] <- "HAVCR1"
  LR[LR == "TMEM8A"] <- "PGAP6"
  LR[LR == "YARS"] <- "YARS1"
  setnames(
    x = LR,
    old = c("ligand", "receptor"),
    new = c("L1", "R1")
  )
  #remove genes that should be excluded according to our manual curation
  LR <- LR[
    !(L1 %in% GENES_to_remove_human$gene) &
      !(R1 %in% GENES_to_remove_human$gene)
  ]
  if (species %in% c("mouse", "rat")) {
    ortho <- get_orthologs(
      genes = unique(c(LR$L1, LR$R1)),
      input_species = "human",
      output_species = species,
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    LR <- merge_LR_orthologs(
      LR_dt = LR,
      ortho_dt = ortho,
      nL = 1,
      charL = "L",
      nR = 1,
      charR = "R",
      input_species = "human",
      output_species = species
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
  LR[, source := gsub("fantom5", "FANTOM5", source)]
  LR[, source := gsub("literature;", "", source)]
  LR[, source := gsub("literature", "", source)]
  LR[, source := gsub("uniprot", "PPI", source)]
  LR[, source := ifelse(source == "", NA, source)]
  LR[, SOURCE := paste(PMIDs, source, sep = ";")]
  LR[, SOURCE := gsub(";NA", "", SOURCE)]
  LR[, SOURCE := gsub("NA;", "", SOURCE)]
  return(
    list(
      retrieved_date = retrieved_date,
      retrieved_from = retrieved_from,
      LR = LR[, cols_to_keep, with = FALSE]
    )
  )
}

prepare_LR_NicheNet <- function(
  species,
  one2one
) {
  LR_SORTED <- SOURCE <- L1 <- R1 <- NULL
  # if (!requireNamespace("nichenetr", quietly = TRUE)) {
  #   stop(
  #     paste0(
  #       "Package \"nichenetr\" needed for this function to work.",
  #       "Please install it."
  #     ),
  #     call. = FALSE
  #   )
  # }
  retrieved_date <- as.Date("2021-03-22")
  retrieved_from <- "nichenetr::lr_network"
  #LR <- nichenetr::lr_network
  LR <- NICHENET_human
  setDT(LR)
  #there are 8 genes that are not HGNC approved symbol, we change them
  LR[LR == "ATP5B"] <- "ATP5F1B"
  LR[LR == "CTGF"] <- "CCN2"
  LR[LR == "CYR61"] <- "CCN1"
  LR[LR == "HFE2"] <- "HJV"
  LR[LR == "NOV"] <- "CCN3"
  LR[LR == "WISP2"] <- "CCN5"
  LR[LR == "WISP3"] <- "CCN6"
  LR[LR == "YARS"] <- "YARS1"
  setnames(
    x = LR,
    old = c("from", "to", "source"),
    new = c("L1", "R1", "SOURCE")
  )
  #remove genes that should be excluded according to our manual curation
  LR <- LR[
    !(L1 %in% GENES_to_remove_human$gene) &
      !(R1 %in% GENES_to_remove_human$gene)
  ]
  if (species %in% c("mouse", "rat")) {
    ortho <- get_orthologs(
      genes = unique(c(LR$L1, LR$R1)),
      input_species = "human",
      output_species = species,
      one2one = one2one
    )
    ortho <- stats::na.omit(ortho)
    LR <- merge_LR_orthologs(
      LR_dt = LR,
      ortho_dt = ortho,
      nL = 1,
      charL = "L",
      nR = 1,
      charR = "R",
      input_species = "human",
      output_species = species
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
  source_temp_kegg <- paste0(
    source_temp[grepl("kegg", source_temp)],
    collapse = "|"
  )
  source_temp_ppi <- paste0(
    source_temp[grepl("ppi", source_temp)],
    collapse = "|"
  )
  LR[, SOURCE := gsub(source_temp_ppi, "PPI", SOURCE)]
  LR[, SOURCE := gsub(source_temp_kegg, "KEGG:NicheNet", SOURCE)]
  LR[, SOURCE := gsub("pharmacology", "IUPHAR", SOURCE)]
  LR[, SOURCE := gsub("ramilowski_known", "FANTOM5", SOURCE)]
  return(
    list(
      retrieved_date = retrieved_date,
      retrieved_from = retrieved_from,
      LR = LR[, cols_to_keep, with = FALSE]
    )
  )
}

# prepare_LR_scTensor <- function(
#   species
# ) {
#   SOURCE <- NULL
#   if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
#     stop(
#       paste0(
#         "Package \"AnnotationDbi\" needed for this function to work.",
#         "Please install it."
#       ),
#       call. = FALSE
#     )
#   }
#   if (!requireNamespace("LRBase.Mmu.eg.db", quietly = TRUE)) {
#     stop(
#       paste0(
#
#         "Package \"LRBase.Mmu.eg.db\" needed for this function to work.",
#         "Please install it."
#       ),
#       call. = FALSE
#     )
#   }
#   if (!requireNamespace("LRBase.Hsa.eg.db", quietly = TRUE)) {
#     stop(
#       paste0(
#         "Package \"LRBase.Hsa.eg.db\" needed for this function to work.",
#         "Please install it."
#       ),
#       call. = FALSE
#     )
#   }
#   if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
#     stop(
#       paste0(
#         "Package \"AnnotationHub\" needed for this function to work.",
#         "Please install it."
#       ),
#       call. = FALSE
#     )
#   }
#   LR_SORTED <- LIGAND_1 <- RECEPTOR_1 <- NULL
#   retrieved_date <- as.Date("2021-03-22")
#   retrieved_from <- "LRBaseDbi"
#   if (species %in% c("mouse", "rat")) {
#     key <- AnnotationDbi::keys(
#       LRBase.Mmu.eg.db::LRBase.Mmu.eg.db,
#       keytype = "GENEID_L"
#     )
#     LR <- AnnotationDbi::select(
#       LRBase.Mmu.eg.db::LRBase.Mmu.eg.db,
#       keys = key,
#       columns = c("GENEID_L", "GENEID_R", "SOURCEDB"),
#       keytype = "GENEID_L"
#     )
#   }
#   if (species == "human"){
#     key <- AnnotationDbi::keys(
#       LRBase.Hsa.eg.db::LRBase.Hsa.eg.db,
#       keytype = "GENEID_L"
#     )
#     LR <- AnnotationDbi::select(
#       LRBase.Hsa.eg.db::LRBase.Hsa.eg.db,
#       keys = key,
#       columns = c("GENEID_L", "GENEID_R", "SOURCEDB"),
#       keytype = "GENEID_L"
#     )
#   }
#   LR <- unique(LR)
#   ah <- AnnotationHub::AnnotationHub()
#   if (species == "mouse") {
#     hs <- AnnotationHub::query(
#       ah,
#       c("OrgDb", "Mus musculus")
#     )[[1]]
#   }
#   if (species == "human") {
#     hs <- AnnotationHub::query(
#       ah,
#       c("OrgDb", "Homo sapiens")
#     )[[1]]
#   }
#   LR_L_match <- AnnotationDbi::select(
#     hs,
#     column = c("SYMBOL", "ENTREZID"),
#     keytype = "ENTREZID",
#     keys = as.character(LR$GENEID_L)
#   )
#   LR_R_match <- AnnotationDbi::select(
#     hs,
#     column = c("SYMBOL", "ENTREZID"),
#     keytype = "ENTREZID",
#     keys = as.character(LR$GENEID_R)
#   )
#   if (identical(LR_L_match$ENTREZID, as.character(LR$GENEID_L))) {
#     LR$LIGAND_1 <- LR_L_match$SYMBOL
#   } else {
#     stop("Matching not possible.")
#   }
#   if (identical(LR_R_match$ENTREZID, as.character(LR$GENEID_R))) {
#     LR$RECEPTOR_1 <- LR_R_match$SYMBOL
#   } else {
#     stop("Matching not possible.")
#   }
#   setDT(LR)
#   LR[, c("GENEID_L", "GENEID_R") := list(NULL, NULL)]
#   setnames(LR, old = "SOURCEDB", new = "SOURCE")
#   LR[, LR_SORTED := list(sapply(1:nrow(.SD), function(i) {
#     temp <- c(LIGAND_1[[i]], RECEPTOR_1[[i]])
#     temp <- temp[!is.na(temp)]
#     temp <- sort(temp)
#     temp <- paste0(temp, collapse = "_")
#   }))]
#   LR <- LR[!duplicated(LR_SORTED)]
#   LR <- LR[!is.na(LIGAND_1) & !is.na(RECEPTOR_1)]
#   cols_to_keep <- c(
#     "LR_SORTED", "SOURCE",
#     "LIGAND_1", "RECEPTOR_1"
#   )
#   #remove genes that should be excluded according to our manual curation
#   if (species == "mouse") {
#     LR <- LR[
#       !(LIGAND_1 %in% GENES_to_remove_mouse$gene) &
#         !(RECEPTOR_1 %in% GENES_to_remove_mouse$gene)
#     ]
#   }
#   if (species == "human") {
#     LR <- LR[
#       !(LIGAND_1 %in% GENES_to_remove_human$gene) &
#         !(RECEPTOR_1 %in% GENES_to_remove_human$gene)
#     ]
#   }
#   if (species == "mouse") {
#     LR[, SOURCE := gsub("ENSEMBL_DLRP|NCBI_DLRP", "SCT:DLRP", SOURCE)]
#     LR[, SOURCE := gsub("ENSEMBL_IUPHAR|NCBI_IUPHAR", "SCT:IUPHAR", SOURCE)]
#     LR[, SOURCE := gsub("ENSEMBL_HPMR|NCBI_HPMR", "SCT:HPMR", SOURCE)]
#     LR[, SOURCE := gsub(
#       "ENSEMBL_CELLPHONEDB|NCBI_CELLPHONEDB",
#       "SCT:CellPhoneDB",
#       SOURCE
#     )]
#     LR[, SOURCE := gsub(
#       "ENSEMBL_SINGLECELLSIGNALR|NCBI_SINGLECELLSIGNALR",
#       "SCT:SingleCellSignalR",
#       SOURCE
#     )]
#     LR[, SOURCE := gsub("SWISSPROT_STRING", "PPI", SOURCE)]
#     LR[, SOURCE := gsub("TREMBL_STRING", "PPI", SOURCE)]
#   }
#   if (species == "human") {
#     LR[, SOURCE := gsub(
#       "SWISSPROT_STRING|TREMBL_STRING|BADERLAB",
#       "PPI",
#       SOURCE)]
#     LR[, SOURCE := gsub("SWISSPROT_HPRD|TREMBL_HPRD", "HPRD", SOURCE)]
#     LR[, SOURCE := gsub("IUPHAR", "SCT:IUPHAR", SOURCE)]
#     LR[, SOURCE := gsub("DLRP", "SCT:DLRP", SOURCE)]
#     LR[, SOURCE := gsub("HPMR", "SCT:HPMR", SOURCE)]
#     LR[, SOURCE := gsub("FANTOM5", "SCT:FANTOM5", SOURCE)]
#     LR[, SOURCE := gsub("CELLPHONEDB", "SCT:CellPhoneDB", SOURCE)]
#     LR[, SOURCE := gsub("SINGLECELLSIGNALR", "SCT:SingleCellSignalR", SOURCE)]
#   }
#   return(
#     list(
#       retrieved_date = retrieved_date,
#       retrieved_from = retrieved_from,
#       LR = LR[, cols_to_keep, with = FALSE]
#     )
#   )
# }

get_KEGG_PWS_interactions <- function(
  species,
  LR_db
) {
  KEGG_NAME <- GENE <- i.KEGG_NAME <- NULL
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"KEGGREST\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (species == "mouse") {
    organism <- "mmu"
    string_to_remove <- " - Mus musculus (mouse)"
  }
  if (species == "human") {
    organism <- "hsa"
    string_to_remove <-  " - Homo sapiens (human)"
  }
  if (species == "rat") {
    organism <- "rno"
    string_to_remove <-  " - Rattus norvegicus (rat)"
  }
  KEGG_PW <- KEGGREST::keggList(
    database = "pathway",
    organism = organism
  )
  KEGG_PW <- data.table(
    "KEGG_ID" = names(KEGG_PW),
    "KEGG_NAME" = KEGG_PW
  )
  KEGG_PW[, KEGG_NAME := gsub(string_to_remove, "", KEGG_NAME, fixed = TRUE)]
  KEGG_PW_to_genes <- rbindlist(
    l = lapply(
      KEGG_PW$KEGG_ID,
      function(id) {
        temp <- KEGGREST::keggGet(
          dbentries = id
        )
        temp <- temp[[1]]$GENE
        temp <- temp[grepl(";", temp)]
        temp <- sort(gsub(";.*", "", temp))
        if (length(temp) > 0) {
          result <- data.table(
            GENE = temp,
            KEGG_ID = id
          )
        } else {
          result <- NULL
        }
        return(result)
      }
    )
  )
  LR_KEGG_PW <- rbindlist(
    apply(
      LR_db,
      MARGIN = 1,
      function(row) {
        LIGAND_PW <- unique(KEGG_PW_to_genes[
          GENE %in% c(row[["LIGAND_1"]], row[["LIGAND_2"]])
        ]$KEGG_ID)
        RECEPTOR_PW <- unique(KEGG_PW_to_genes[
          GENE %in% c(
            row[["RECEPTOR_1"]],
            row[["RECEPTOR_2"]],
            row[["RECEPTOR_3"]]
          )
        ]$KEGG_ID)
        res_inter <- intersect(LIGAND_PW, RECEPTOR_PW)
        if (length(res_inter) > 0) {
          res_inter <- data.table(
            LRI = rep(row[["LRI"]], length(res_inter)),
            KEGG_ID = res_inter
          )
        } else {
          res_inter <- NULL
        }
        return(res_inter)
      }
    )
  )
  LR_KEGG_PW[
    KEGG_PW,
    on = "KEGG_ID",
    KEGG_NAME := i.KEGG_NAME
  ]
  return(LR_KEGG_PW)
}

get_GO_interactions <- function(
  species,
  LR_db,
  only_genes_annotations = FALSE
) {
  GO_NAME <- i.GO_name <- ASPECT <- GO_ID <-
    LEVEL <- i.LEVEL <- NULL
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"biomaRt\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("ontoProc", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"ontoProc\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("ontologyIndex", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"ontologyIndex\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  ALL_LR_genes <- unique(
    unlist(
      LR_db[
        ,
        c("LIGAND_1", "LIGAND_2", "RECEPTOR_1",
          "RECEPTOR_2", "RECEPTOR_3")
      ]
    )
  )
  ALL_LR_genes <- ALL_LR_genes[!is.na(ALL_LR_genes)]
  if (species == "mouse") {
    dataset <- "mmusculus_gene_ensembl"
    id_gene <- "mgi_symbol"
  }
  if (species == "human") {
    dataset <- "hsapiens_gene_ensembl"
    id_gene <- "hgnc_symbol"
  }
  if (species == "rat") {
    dataset <- "rnorvegicus_gene_ensembl"
    id_gene <- "rgd_symbol"
  }
  mart <- biomaRt::useMart(
    "ensembl",
    host = "https://nov2020.archive.ensembl.org",
    dataset = dataset,
    verbose = TRUE
  )
  ALL_LR_genes_info <- biomaRt::getBM(
    attributes = c(
      id_gene,
      "go_id",
      "name_1006"
    ),
    filters = id_gene,
    mart = mart,
    values = ALL_LR_genes
  )
  setDT(ALL_LR_genes_info)
  if (only_genes_annotations) {
    return(ALL_LR_genes_info)
  }
  onto_go_terms <- ontoProc::getGeneOnto()
  go_names <- onto_go_terms$name
  ALL_LR_genes_go <- sapply(
    ALL_LR_genes,
    function(gene) {
      temp_go <- unique(ALL_LR_genes_info[get(id_gene) == gene]$name_1006)
      ontologyIndex::get_ancestors(
        onto_go_terms,
        names(go_names[go_names %in% temp_go])
      )
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  )
  LR_interactions_go_intersection <- rbindlist(
    apply(
      LR_db,
      MARGIN = 1,
      function(row) {
        LIGAND_GO <- unique(c(
          ALL_LR_genes_go[[row[["LIGAND_1"]]]],
          ALL_LR_genes_go[[row[["LIGAND_2"]]]]
        ))
        RECEPTOR_GO <- unique(c(
          ALL_LR_genes_go[[row[["RECEPTOR_1"]]]],
          ALL_LR_genes_go[[row[["RECEPTOR_2"]]]],
          ALL_LR_genes_go[[row[["RECEPTOR_3"]]]]
        ))
        res_inter <- intersect(LIGAND_GO, RECEPTOR_GO)
        if (length(res_inter) > 0) {
          res_inter <- data.table(
            LRI = rep(row[["LRI"]], length(res_inter)),
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
  LR_interactions_go_intersection[
    go_id_name_dt,
    on = "GO_ID==ID",
    GO_NAME := i.GO_name
  ]
  return(LR_interactions_go_intersection)
}

get_ECM_genes <- function(
  species
) {
  GO_ID <- mgi_symbol <- NULL
  LRI_curated = scDiffCom::LRI_mouse$LRI_curated
  GO_interactions = get_GO_interactions(
    species,
    LRI_curated,
    only_genes_annotations = TRUE
  )
  GO_interactions = GO_interactions[GO_ID != ""]

  get_ECM_GOs <- function() {
    return(
      c(
        "GO:0031012", #
        "GO:0005578",
        "GO:0005201", #
        "GO:1990430",
        "GO:0035426" # Found in LRI_mouse$go curated,not in get_GO_interactions
      )
    )
  }
  ECM_GOs = get_ECM_GOs()
  ecm_genes = GO_interactions[GO_ID %in% ECM_GOs, mgi_symbol]
  return(ecm_genes)
}

merge_LR_orthologs <- function(
  LR_dt,
  ortho_dt,
  nL,
  charL,
  nR,
  charR,
  input_species,
  output_species
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
  merge_id <- c(paste0(output_species, "_symbol"), "confidence", "type")
  LR_temp <- copy(LR_dt)
  LR_temp[
    ,
    c(out_names) :=
      c(
        sapply(
          1:nL,
          function(i) {
            as.list(
              ortho_dt[
                .SD,
                on = c(paste0(input_species, "_symbol==", charL, i)),
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
                       on = c(paste0(input_species, "_symbol==", charR, i)),
                       mget(paste0("x.", merge_id))
              ]
            )
          }
        )
      )
  ]
  LR_temp <- stats::na.omit(LR_temp, cols = c("LIGAND_1", "RECEPTOR_1"))
  LR_temp[, to_keep := sapply(1:nrow(.SD), function(i) {
    all(c(
      sapply(1:nL, function(j) {
        !(is.na(get(paste0("LIGAND_", j))[[i]]) &
            !is.na(get(paste0(charL, j))[[i]]))
      }),
      sapply(1:nR, function(j) {
        !(is.na(get(paste0("RECEPTOR_", j))[[i]]) &
            !is.na(get(paste0(charR, j))[[i]]))
      })
    ))
  })]
  LR_temp <- LR_temp[to_keep == TRUE]
  LR_temp[, to_keep := NULL]
  if (output_species == "mouse") {
    #there are genes that are not MGI approved symbol, we change them
    LR_temp[LR_temp == "Oit1"] <- "Fam3d"
    LR_temp[LR_temp == "Il1f5"] <- "Il36rn"
    LR_temp[LR_temp == "Il1f6"] <- "Il36a"
    LR_temp[LR_temp == "Il1f8"] <- "Il36b"
  }
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
  output_species,
  one2one
) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"biomaRt\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  ensembl_gene_id <- inl <- outl <- output <- input <- confidence <- NULL
  if (input_species == "mouse") {
    id_in <- "mmusculus"
    id_gene <- "mgi_symbol"
    name_in <- "mouse"
  } else if (input_species == "human") {
    id_in <- "hsapiens"
    id_gene <- "hgnc_symbol"
    name_in <- "human"
  } else if (input_species == "rat") {
    id_in <- "rnorvegicus"
    id_gene <- "rgd_symbol"
    name_in <- "rat"
  } else {
    stop("Input species not supported in function get_orthologs.")
  }
  if (output_species == "mouse") {
    id_out <- "mmusculus"
    name_out <- "mouse"
  } else if (output_species == "human") {
    id_out <- "hsapiens"
    name_out <- "human"
  } else if (output_species == "rat") {
    id_out <- "rnorvegicus"
    name_out <- "rat"
  } else {
    stop("Output species not supported in function get_orthologs.")
  }
  dataset <- paste0(id_in, "_gene_ensembl")
  gene_name <- paste0(id_out, "_homolog_associated_gene_name")
  ortho_confidence <- paste0(id_out, "_homolog_orthology_confidence")
  ortho_type <- paste0(id_out, "_homolog_orthology_type")
  mart <- biomaRt::useMart(
    "ensembl",
    host = "https://nov2020.archive.ensembl.org",
    dataset = dataset,
    verbose = TRUE
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
    ensembl_all <- ensembl_all[
      eval(as.symbol(ortho_type)) == "ortholog_one2one",
    ]
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
      warning(
        paste0(
          "There are some duplicates from orthology conversion.",
          "Removing them by using 'unique'."
        )
      )
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

get_GO_LEVELS <- function(
) {
  LEVEL <- N_ANCESTORS <- ID <- ASPECT <- NULL
  if (!requireNamespace("ontoProc", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"ontoProc\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("ontologyIndex", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"ontologyIndex\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  ontoGO <- ontoProc::getGeneOnto()
  GO_summary <- data.table(
    ID = ontoGO$id,
    NAME = ontoGO$name,
    N_ANCESTORS = sapply(
      ontoGO$ancestors,
      length
    ),
    N_PARENTS = sapply(
      ontoGO$parents,
      length
    ),
    N_CHILDREN = sapply(
      ontoGO$children,
      length
    )
  )
  GO_summary[, LEVEL := ifelse(
    N_ANCESTORS == 1,
    1,
    Inf
  )]
  GO_summary[, LEVEL := ifelse(
    N_ANCESTORS == 2,
    2,
    LEVEL
  )]
  find_level_n <- function(
    n,
    go_table
  ) {
    LEVEL <- N_CHILDREN <- NULL
    if (nrow(go_table[LEVEL == n - 1]) == 0) {
      stop("Level n-1 must be performed before level n")
    }
    temp_n <- unique(
      unlist(
        ontoGO$children[
          go_table[
            LEVEL == n - 1 & N_CHILDREN > 0
          ]$ID
        ]
      )
    )
    temp_n_ancestors <- ontoGO$ancestors[temp_n]
    temp_n_logical <- sapply(
      seq_along(temp_n_ancestors),
      function(i) {
        all(
          temp_n_ancestors[[i]] %in% c(
            go_table[LEVEL <= n - 1]$ID,
            names(temp_n_ancestors)[[i]]
          )
        )
      }
    )
    names(temp_n_logical) <- names(temp_n_ancestors)
    return(names(temp_n_logical[temp_n_logical]))
  }
  i <- 3
  while (Inf %in% GO_summary$LEVEL) {
    message(paste0("Processing level ", i))
    GO_summary[, LEVEL := ifelse(
      ID %in% find_level_n(i, .SD),
      i,
      LEVEL
    )]
    i <- i + 1
  }
  all_bp_desc <- ontologyIndex::get_descendants(ontoGO, "GO:0008150")
  all_mf_desc <- ontologyIndex::get_descendants(ontoGO, "GO:0003674")
  all_cc_desc <- ontologyIndex::get_descendants(ontoGO, "GO:0005575")
  GO_summary[, ASPECT := ifelse(
    ID %in% all_bp_desc,
    "biological_process",
    ifelse(
      ID %in% all_mf_desc,
      "molecular_function",
      ifelse(
        ID %in% all_cc_desc,
        "cellular_component",
        "other"
      )
    )
  )]
  return(GO_summary)
}
