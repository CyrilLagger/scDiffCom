build_interactive_network <- function(
  object,
  network_type,
  layout_type,
  class_signature,
  subobject_name,
  abbreviation_table
) {
  ID <- VALUE <- i.ABBR_CELLTYPE <- EMITTER_CELLTYPE <-
    RECEIVER_CELLTYPE <- NULL
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"visNetwork\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"RColorBrewer\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"kableExtra\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"igraph\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop(
      paste0(
        "Package \"igraph\" needed for this function to work.",
        "Please install it."
      ),
      call. = FALSE
    )
  }
  if (!object@parameters$permutation_analysis) {
    stop(
      paste0(
        "No network available for the current object: ",
        "permutation analysis must have been performed"
      )
    )
  }
  if (!object@parameters$conditional_analysis &
      network_type != "condition1_network") {
    stop(
      paste0(
        "No network of type ",
        network_type,
        " available for the current object: ",
        "differential analysis must have been performed",
        " or select `condition1_network` as network_type"
      )
    )
  }
  cci_table_detected <- copy(object@cci_table_detected)
  if (class_signature == "scDiffComCombined") {
    object_name <- subobject_name
    cci_table_detected <- cci_table_detected[
      ID == subobject_name
      ][, ID := NULL]
  } else {
    object_name <- object@parameters$object_name
  }
  if(!is.null(abbreviation_table)) {
    if(!methods::is(abbreviation_table, "data.table") ||
       !methods::is(abbreviation_table, "data.frame") ||
       !identical(
         names(abbreviation_table),
         c("ORIGINAL_CELLTYPE", "ABBR_CELLTYPE")
       )) {
      warning(
        paste0(
          "No abbreviation will be used:",
          " `abbreviation table` must be a 2-colum data.frame or data.table",
          "with names ORIGINAL_CELLTYPE and ABBR_CELLTYPE"
          )
        )
      abbreviation_table <- NULL
    } else {
      setDT(abbreviation_table)
      actual_celltypes <- union(
        cci_table_detected[["EMITTER_CELLTYPE"]],
        cci_table_detected[["RECEIVER_CELLTYPE"]]
      )
      if (!identical(
        sort(actual_celltypes),
        sort(abbreviation_table[["ORIGINAL_CELLTYPE"]])
      )) {
        warning(
          paste0(
            "No abbreviation will be used:",
            " `abbreviation table` must contain",
            " a column with the original cell-types")
        )
        abbreviation_table <- NULL
      } else if (sum(duplicated(abbreviation_table)) > 0) {
        warning(
          paste0(
            "No abbreviation will be used:",
            " `abbreviation table` must not contain duplicated rows"))
        abbreviation_table <- NULL
      } else {
        cci_table_detected[
          ,
          "EMITTER_CELLTYPE_ORIGINAL" := EMITTER_CELLTYPE
        ]
        cci_table_detected[
          ,
          "RECEIVER_CELLTYPE_ORIGINAL" := RECEIVER_CELLTYPE
        ]
        cci_table_detected[
          abbreviation_table,
          on = "EMITTER_CELLTYPE==ORIGINAL_CELLTYPE",
          "EMITTER_CELLTYPE" := i.ABBR_CELLTYPE]
        cci_table_detected[
          abbreviation_table,
          on = "RECEIVER_CELLTYPE==ORIGINAL_CELLTYPE",
          "RECEIVER_CELLTYPE" := i.ABBR_CELLTYPE]
      }
    }
  }
  if (!object@parameters$conditional_analysis) {
    network <- interactive_from_igraph(
      cci_table_detected = cci_table_detected,
      conds = NULL,
      ora_table_ER = NULL,
      ora_table_LR = NULL,
      network_type = "condition1_network",
      layout_type = layout_type,
      object_name = object_name
    )
  } else {
    conds <- c(
      object@parameters$seurat_condition_id$cond1_name,
      object@parameters$seurat_condition_id$cond2_name
    )
    if (network_type == "ORA_network") {
      if (identical(object@ora_table, list()) |
          is.null(object@ora_table$ER_CELLTYPES) |
          is.null(object@ora_table$LRI)) {
        stop(paste0(
          "No network of type ",
          network_type,
          " available for the current object: ",
          " the object must contain an ORA table"
        ))
      }
      ora_table_ER <- copy(object@ora_table$ER_CELLTYPES)
      ora_table_LR <- copy(object@ora_table$LRI)
      if (class_signature == "scDiffComCombined") {
        ora_table_ER <- ora_table_ER[ID == subobject_name][, ID := NULL]
        ora_table_LR <- ora_table_LR[ID == subobject_name][, ID := NULL]
      }
      ora_table_ER[, c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE") := list(
        sub("_.*", "", VALUE),
        sub(".*_", "", VALUE)
      )]
      if (!is.null(abbreviation_table)) {
        ora_table_ER[
          abbreviation_table,
          on = "EMITTER_CELLTYPE==ORIGINAL_CELLTYPE",
          "EMITTER_CELLTYPE" := i.ABBR_CELLTYPE]
        ora_table_ER[
          abbreviation_table,
          on = "RECEIVER_CELLTYPE==ORIGINAL_CELLTYPE",
          "RECEIVER_CELLTYPE" := i.ABBR_CELLTYPE]
      }
    } else {
      ora_table_ER <- NULL
      ora_table_LR <- NULL
    }
    network <- interactive_from_igraph(
      cci_table_detected = cci_table_detected,
      conds = conds,
      ora_table_ER = ora_table_ER,
      ora_table_LR = ora_table_LR,
      network_type = network_type,
      layout_type = layout_type,
      object_name = object_name
    )
  }
  return(network)
}

interactive_from_igraph <- function(
  cci_table_detected,
  conds,
  ora_table_ER,
  ora_table_LR,
  network_type,
  layout_type,
  object_name
) {
  config <- setup_graph_config()
  G <- build_igraph(
    cci_table_detected = cci_table_detected,
    conds = conds,
    ora_table_ER = ora_table_ER,
    ora_table_LR = ora_table_LR,
    network_type = network_type,
    layout_type = layout_type,
    config
  )
  interactive_network <- build_visnetwork(
    G = G,
    object_name = object_name,
    network_type = network_type,
    layout_type = layout_type,
    config = config
  )
  return(interactive_network)
}

setup_graph_config <- function(
) {
  GRAPH_CONFIG <- list(
    EDGE_COLORING = list(
      ORA_COLOR_UP = "#F94144", # red
      ORA_COLOR_DOWN = "#277DA1", # blue
      ORA_COLOR_DIFF = "#F9C74F",
      ORA_COLOR_FLAT = "#90BE6D",
      ORA_COLOR_NONE = grDevices::rgb(0.2, 0.2, 0.2, alpha = 0.1),
      BREWER_N = 7,
      BREWER_NAME = "RdBu"
    ),
    ORA_PARAMETERS = list(
      CUTOFF_OR = 1,
      CUTOFF_BHP = 0.05
    ),
    NODE_COLORING = list(
      BACKGROUND = "#f8961e",
      BORDER = "#577590",
      HIGHLIGHT = list(
        BACKGROUND = "#43AA8B",
        BORDER = "#577590"
      ),
      HOVER = list(
        BACKGROUND = "#43AA8B",
        BORDER = "#577590"
      )
    ),
    EDGE_STYLE = list(
      WIDTH_MIN = 0,
      WIDTH_MAX = 8,
      WIDTH = 2.5,
      ARROW_SIZE = 0.5
    ),
    LAYOUT = list(
      HGAP = 20,
      VGAP = 20,
      DISPERSE = FALSE,
      IGRAPH_FUN = "circle" #nicely, sugiyama, with_fr
    ),
    VERTEX_STYLE = list(
      MAXSIZE = 5, #20,
      MINSIZE = 5, #5,
      SIZE = 5,
      LABEL_DIST = 1.5, #1.5,
      LABEL_CEX = 1.5, # 1.2,
      COLOR = "#33FF66"
    ),
    LEGEND = list(
      LEGEND_LABELS = c(
        "Significant, but small effect",
        "Upregulated",
        "Downregulated",
        "Altered"
      ),
      PCH = c(15),
      CEX = 0.7,
      PT.CEX = 1,
      BG = "#CCCCCC",
      NCOL = 2
    ),
    VISNETWORK = list(
      WIDTH = "100%",
      HEIGHT = NULL,
      BACKGROUND = 	"#F5F5F5"
    )
  )
  return(GRAPH_CONFIG)
}

build_igraph <- function(
  cci_table_detected,
  conds,
  ora_table_ER,
  ora_table_LR,
  network_type,
  layout_type,
  config = config
) {
  edge_table <- build_edge_table(
    cci_table_detected = cci_table_detected,
    conds = conds,
    ora_table_ER = ora_table_ER,
    ora_table_LR = ora_table_LR,
    network_type = network_type,
    layout_type = layout_type,
    config = config
  )
  vertex_table <- build_vertex_table(
    cci_table_detected = cci_table_detected,
    edge_table = edge_table,
    conds = conds,
    network_type = network_type,
    layout_type = layout_type,
    config = config
  )
  G <- igraph::graph_from_data_frame(
    d = edge_table,
    directed = TRUE,
    vertices = vertex_table
  )
  G <- setup_layout(
    G = G,
    network_type,
    layout_type,
    config = config
  )
  return(G)
}

build_edge_table <- function(
  cci_table_detected,
  conds,
  ora_table_ER,
  ora_table_LR,
  network_type,
  layout_type,
  config
) {
  ORA_TYPE <- NULL
  edge_table <- extract_edge_metadata(
    cci_table_detected = cci_table_detected,
    conds = conds,
    ora_table_ER = ora_table_ER,
    network_type = network_type,
    config = config
  )
  edge_table <- add_edge_layout(
    edge_table = edge_table,
    conds = conds,
    ora_table_LR = ora_table_LR,
    network_type = network_type,
    layout_type,
    config = config
  )
  if (network_type == "ORA_network") {
    edge_table <- edge_table[ORA_TYPE != "NONE"]
  }
  return(edge_table)
}

extract_edge_metadata <- function(
  cci_table_detected,
  conds,
  ora_table_ER,
  network_type,
  config
) {
  IS_CCI_EXPRESSED <- i.N <- REGULATION <-
    NUM_LRIS_UP <- NUM_LRIS_DOWN <- i.ORA_TYPE <- NULL
  all_cell_types <- union(
    unique(cci_table_detected[["EMITTER_CELLTYPE"]]),
    unique(cci_table_detected[["RECEIVER_CELLTYPE"]])
  )
  edge_table <- CJ(
    "from" = all_cell_types,
    "to" = all_cell_types
  )
  if (is.null(conds)) {
    edge_table[
      cci_table_detected[
        IS_CCI_EXPRESSED == TRUE,
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      "NUM_LRIS" := i.N
    ]
    edge_table[is.na(edge_table)] <- 0
  } else {
    edge_table[
      cci_table_detected[
        get(paste0("IS_CCI_EXPRESSED_", conds[[1]])) == TRUE |
          get(paste0("IS_CCI_EXPRESSED_", conds[[2]])) == TRUE,
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      "NUM_LRIS_TOTAL" := i.N
    ]
    edge_table[
      cci_table_detected[
        get(paste0("IS_CCI_EXPRESSED_", conds[[1]])) == TRUE,
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      paste0("NUM_LRIS_", conds[[1]]) := i.N
    ]
    edge_table[
      cci_table_detected[
        get(paste0("IS_CCI_EXPRESSED_", conds[[2]])) == TRUE,
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      paste0("NUM_LRIS_", conds[[2]]) := i.N
    ]
    edge_table[is.na(edge_table)] <- 0
    edge_table[, "NUM_LRIS_REL_DIFF" :=
                 get(paste0("NUM_LRIS_", conds[[2]])) -
                 get(paste0("NUM_LRIS_", conds[[1]]))]
    edge_table[
      cci_table_detected[
        REGULATION == "UP",
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      "NUM_LRIS_UP" := i.N
    ]
    edge_table[
      cci_table_detected[
        REGULATION == "DOWN",
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      "NUM_LRIS_DOWN" := i.N
    ]
    edge_table[
      cci_table_detected[
        REGULATION == "FLAT",
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      "NUM_LRIS_FLAT" := i.N
    ]
    edge_table[
      cci_table_detected[
        REGULATION == "NSC",
        .N,
        by = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE")],
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      "NUM_LRIS_NSC" := i.N
    ]
    edge_table[, "NUM_LRIS_DIFF" := NUM_LRIS_UP + NUM_LRIS_DOWN]
  }
  if (network_type == "ORA_network") {
    cols_to_add <- c(
      "ORA_TYPE",
      "OR_UP", "OR_DOWN", "OR_FLAT",
      "BH_P_VALUE_UP", "BH_P_VALUE_DOWN", "BH_P_VALUE_FLAT"
    )
    edge_table[
      process_celltype_pairs_enrichment(
        ora_table_ER = ora_table_ER,
        config = config
      ),
      on = c("from==EMITTER_CELLTYPE", "to==RECEIVER_CELLTYPE"),
      (cols_to_add) := mget(paste0("i.", cols_to_add))
    ]
  }
  # TODO: Enrich with other info: cell families, GOs, ORA on cell types.
  return(edge_table)
}

process_celltype_pairs_enrichment <- function(
  ora_table_ER,
  config
) {
  OR_UP <- OR_DOWN <- OR_FLAT <-
    BH_P_VALUE_UP <- BH_P_VALUE_DOWN <- BH_P_VALUE_FLAT <- NULL
  dt_ora <- copy(ora_table_ER)
  OR_MIN <- config$ORA_PARAMETERS$CUTOFF_OR
  BH_MAX <- config$ORA_PARAMETERS$CUTOFF_BHP
  dt_ora[, "ORA_TYPE" := ifelse(
    OR_UP >= OR_MIN &
      BH_P_VALUE_UP <= BH_MAX &
      OR_DOWN >= OR_MIN &
      BH_P_VALUE_DOWN <= BH_MAX,
    "DIFF",
    ifelse(
      OR_UP >= OR_MIN & BH_P_VALUE_UP <= BH_MAX,
      "UP",
      ifelse(
        OR_DOWN >= OR_MIN & BH_P_VALUE_DOWN <= BH_MAX,
        "DOWN",
        ifelse(
          OR_FLAT >= OR_MIN & BH_P_VALUE_FLAT <= BH_MAX,
          "FLAT",
          "NONE"
        )
      )
    )
  )]
  cols_to_select1 <- c(
    "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE", "ORA_TYPE"
  )
  cols_to_select2 <- c(
    "OR_UP", "OR_DOWN", "OR_FLAT",
    "BH_P_VALUE_UP", "BH_P_VALUE_DOWN", "BH_P_VALUE_FLAT"
  )
  cols_to_select <- c(cols_to_select1, cols_to_select2)
  dt_ora <- dt_ora[, cols_to_select, with = FALSE]
  dt_ora[
    ,
    (cols_to_select2) := lapply(.SD, signif, 3), .SDcol = cols_to_select2
    ]
  # if (sum(is.na(dt_ora)) > 0 | sum(dt_ora == Inf) > 0 ) {
  #   stop("Inf or NA in `dt_ora`")
  # }
  return(dt_ora)
}

add_edge_layout <- function(
  edge_table,
  conds,
  ora_table_LR,
  network_type,
  layout_type,
  config
) {
  color.color <- i.color <- from <- to <-  NULL
  if (network_type == "ORA_network") {
    edge_table[
      data.table(
        ORA_TYPE = c("DIFF", "UP", "DOWN", "FLAT", "NONE"),
        color = c(
          config$EDGE_COLORING$ORA_COLOR_DIFF,
          config$EDGE_COLORING$ORA_COLOR_UP,
          config$EDGE_COLORING$ORA_COLOR_DOWN,
          config$EDGE_COLORING$ORA_COLOR_FLAT,
          config$EDGE_COLORING$ORA_COLOR_NONE
        )
      ),
      on = "ORA_TYPE" ,
      "color.color" := i.color
    ]
    lab <- "NUM_LRIS_TOTAL"
    main_title <- "Over-represented cell type pairs"
  } else {
    if (network_type == "condition1_network") {
      if (is.null(conds)) {
        lab <- "NUM_LRIS"
        main_title <- paste0(
          "Number of detected ligand-receptor",
          " interactions between cell type pairs")
      } else {
        lab <- paste0("NUM_LRIS_", conds[[1]])
        main_title <- paste0(
          "Number of detected ligand-receptor interactions",
          " between cell type pairs",
          " (", conds[[1]], ")"
        )
      }
    } else if (network_type == "condition2_network") {
      lab <- paste0("NUM_LRIS_", conds[[2]])
      main_title <- paste0(
        "Number of detected ligand-receptor interactions",
        " between cell type pairs",
        " (", conds[[2]], ")"
      )
    } else if (network_type == "up_regulated_network") {
      lab <- "NUM_LRIS_UP"
      main_title <- paste0(
        "Number of up-regulated ligand-receptor interactions",
        " between cell type pairs"
      )
    } else if (network_type == "down_regulated_network") {
      lab <- "NUM_LRIS_DOWN"
      main_title <- paste0(
        "Number of down-regulated ligand-receptor interactions",
        " between cell type pairs"
      )
    } else if (network_type == "difference_network") {
      lab <- "NUM_LRIS_REL_DIFF"
      main_title <- "Difference (Not meaningful!)"
    }
    edge_table[, "color.color" := add_edge_color(get(lab), config)]
    edge_table[, "label" := as.character(get(lab))]
  }
  edge_table[, main_title := main_title]
  edge_table[, "width" := rescale_internal(
    v = sqrt(abs(get(lab))),
    min_ = config$EDGE_STYLE$WIDTH_MIN,
    max_ = config$EDGE_STYLE$WIDTH_MAX
  )]
  edge_table[, "color.highlight" := color.color]
  edge_table[, "color.hover" := color.color]
  edge_table[, "smooth" := TRUE]
  edge_table[, "arrow.size" := config$EDGE_STYLE$ARROW_SIZE]
  edge_table <- edge_annotation_html(edge_table, network_type)
  # num_interacts <- num_interactions_object(
  #   cci_table_detected = cci_table_detected
  #   )
  if (layout_type == "bipartite") {
    edge_table[, "from" := paste0(from, " (E)")]
    edge_table[, "to" := paste0(to, " (R)")]
  }
  if (layout_type == "conventional") {
    edge_table[, "edge.loop.angle" := {
      n <- sqrt(nrow(.SD))
      res <- rep(0, times = n*n)
      temp <- (1:n)*(n+1) - n
      res[seq_along(res) %in% temp] <- rank(-temp)*2*pi/n
      res
    }]
  }
  return(edge_table)
}

add_edge_color <- function(
  NUM_LRIS,
  config
) {
  RColorBrewer::brewer.pal(
    n = config$EDGE_COLORING$BREWER_N,
    name = config$EDGE_COLORING$BREWER_NAME
  )[cut(NUM_LRIS, config$EDGE_COLORING$BREWER_N)]
}

build_vertex_table <- function(
  cci_table_detected,
  edge_table,
  conds,
  network_type,
  layout_type,
  config
) {
  vertex_table <- extract_vertex_metadata(
    cci_table_detected = cci_table_detected,
    conds = conds,
    layout_type = layout_type
  )
  vertex_table <- add_vertex_layout(
    vertex_table = vertex_table,
    edge_table = edge_table,
    conds = conds,
    network_type = network_type,
    layout_type = layout_type,
    config = config
  )
  return(vertex_table)
}

extract_vertex_metadata <- function(
  cci_table_detected,
  conds,
  layout_type
) {
  i.NCELLS_EMITTER <- i.NCELLS_AVG <- name <- vertex_types <- NULL
  all_cell_types <- union(
    unique(cci_table_detected[["EMITTER_CELLTYPE"]]),
    unique(cci_table_detected[["RECEIVER_CELLTYPE"]])
  )
  if ("EMITTER_CELLTYPE_ORIGINAL" %in% colnames(cci_table_detected)) {
    all_cell_types_original <- union(
      unique(cci_table_detected[["EMITTER_CELLTYPE_ORIGINAL"]]),
      unique(cci_table_detected[["RECEIVER_CELLTYPE_ORIGINAL"]])
    )
  }
  if (layout_type == "conventional") {
    if ("EMITTER_CELLTYPE_ORIGINAL" %in% colnames(cci_table_detected)) {
      vertex_table <- data.table(
        name = all_cell_types,
        name_original = all_cell_types_original
      )
    } else {
      vertex_table <- data.table(
        name = all_cell_types
      )
    }
  } else if (layout_type == "bipartite") {
    if ("EMITTER_CELLTYPE_ORIGINAL" %in% colnames(cci_table_detected)) {
      vertex_table <- rbindlist(
        l = list(
          "EMITTER" = data.table(
            name = all_cell_types, #paste0(all_cell_types, " (E)"),
            name_original = all_cell_types_original,
            vertex_types = TRUE
          ),
          "RECEIVER" = data.table(
            name = all_cell_types, # paste0(all_cell_types, " (R)"),
            name_original = all_cell_types_original,
            vertex_types = FALSE
          )
        )
      )
    } else {
      vertex_table <- rbindlist(
        l = list(
          "EMITTER" = data.table(
            name = all_cell_types, #paste0(all_cell_types, " (E)"),
            vertex_types = TRUE
          ),
          "RECEIVER" = data.table(
            name = all_cell_types, # paste0(all_cell_types, " (R)"),
            vertex_types = FALSE
          )
        )
      )
    }
  }
  if (is.null(conds)) {
    NCELLS_TABLE <- unique(cci_table_detected[
      ,
      c("EMITTER_CELLTYPE", "NCELLS_EMITTER"),
      with = FALSE])
    vertex_table[
      NCELLS_TABLE,
      on = "name==EMITTER_CELLTYPE",
      "num_cells" := i.NCELLS_EMITTER
    ]
  } else {
    NCELLS_TABLE <- unique(cci_table_detected[
      ,
      c(
        "EMITTER_CELLTYPE",
        paste0("NCELLS_EMITTER_", conds[[1]]),
        paste0("NCELLS_EMITTER_", conds[[2]])
      ),
      with = FALSE])
    NCELLS_TABLE[, "NCELLS_AVG" :=
                   (get(paste0("NCELLS_EMITTER_", conds[[1]])) +
                      get(paste0("NCELLS_EMITTER_", conds[[2]])))/2]
    vertex_table[
      NCELLS_TABLE,
      on = "name==EMITTER_CELLTYPE",
      "num_cells" := i.NCELLS_AVG
    ]
  }
  if (layout_type == "bipartite") {
    vertex_table[, "name" := ifelse(
      vertex_types,
      paste0(name, " (E)"),
      paste0(name, " (R)")
    )]
  }
  return(vertex_table)
}

add_vertex_layout <- function(
  vertex_table,
  edge_table,
  conds,
  network_type,
  layout_type,
  config
) {
  num_cells <- vertex.size <- name <- vertex_types <- name_original <- NULL
  vertex_table[, "vertex.size" := rescale_internal(
    v = num_cells,
    min_ = config$VERTEX_STYLE$MINSIZE,
    max_ = config$VERTEX_STYLE$MAXSIZE
  )]
  vertex_table[
    ,
    c("id", "label", "value",
      "color.background", "color.border",
      "color.highlight.background", "color.highlight.border",
      "color.hover.background", "color.hover.border",
      "shadow"
    ) := list(
      name, name, vertex.size,
      config$NODE_COLORING$BACKGROUND,
      config$NODE_COLORING$BORDER,
      config$NODE_COLORING$HIGHLIGHT$BACKGROUND,
      config$NODE_COLORING$HIGHLIGHT$BORDER,
      config$NODE_COLORING$HOVER$BACKGROUND,
      config$NODE_COLORING$HOVER$BORDER,
      TRUE
    )
  ]
  if (is.null(conds)) {
    if ("name_original" %in% colnames(vertex_table)) {
      vertex_table[, "title" := paste0(
        "<h3> ", name_original, " </h3><p> Number of Cells: ", num_cells, " </p>"
      )]
    } else {
      vertex_table[, "title" := paste0(
        "<h3> ", name, " </h3><p> Number of Cells: ", num_cells, " </p>"
      )]
    }
  } else {
    if ("name_original" %in% colnames(vertex_table)) {
      vertex_table[, "title" := paste0(
        "<h3> ", name_original, " </h3><p> Avg. Number of Cells: ", num_cells, " </p>"
      )]
    } else {
      vertex_table[, "title" := paste0(
        "<h3> ", name, " </h3><p> Avg. Number of Cells: ", num_cells, " </p>"
      )]
    }
  }
  if (layout_type == "bipartite") {
    vertex_table[, c("group", "level") := list(
      vertex_types,
      ifelse(vertex_types == TRUE, 1, 2)
    )]
    vertex_table[, "vertex_order" := sort_bipartite_vertices(
      vertex_table = .SD,
      edge_table = edge_table,
      network_type = network_type)]
  }
  return(vertex_table)
}

rescale_internal <- function(
  v,
  min_,
  max_
) {
  (v - min(v)) / (max(v) - min(v)) * (max_ - min_) + min_
}

sort_bipartite_vertices <- function(
  vertex_table,
  edge_table,
  network_type
) {
  from <- to <- ORA_TYPE <- vertex_types <- NULL
  if (network_type == "ORA_network") {
    rank_from <- rank(
      sapply(
        vertex_table[vertex_types == TRUE][["name"]],
        function(node) {
          total_edges <- edge_table[from == node, .N]
          total_up <- edge_table[from == node & ORA_TYPE == "UP", .N]
          total_down <- edge_table[from == node & ORA_TYPE == "DOWN", .N]
          total_other <- total_edges - total_down - total_up
          total_down - total_up
        }
      ),
      ties.method = "first"
    )
    rank_to <- rank(
      sapply(
        vertex_table[vertex_types == FALSE][["name"]],
        function(node) {
          total_edges <- edge_table[to == node, .N]
          total_up <- edge_table[to == node & ORA_TYPE == "UP", .N]
          total_down <- edge_table[to == node & ORA_TYPE == "DOWN", .N]
          total_other <- total_edges - total_down - total_up
          total_down - total_up
        }
      ),
      ties.method = "first"
    ) + nrow(vertex_table)/2
  } else {
    # TODO implement special cases for each network_type
    rank_from <- rank(
      (1:nrow(vertex_table[vertex_types == TRUE])),
      ties.method = "first"
    )
    rank_to <- rank(
      (1:nrow(vertex_table[vertex_types == FALSE])),
      ties.method = "first"
    ) +
      nrow(vertex_table)/2
  }
  return(c(rank_from, rank_to))
}

setup_layout <- function(
  G,
  network_type,
  layout_type,
  config
) {
  if (layout_type == "conventional") {
    if (config$LAYOUT$IGRAPH_FUN == "circle") {
      layout <- igraph::layout_in_circle(G)
    } else if (config$LAYOUT$IGRAPH_FUN == "nicely") {
      # determines best layout, likely calls fr
      layout <- igraph::layout_nicely(G, dim=2)
    } else if (config$LAYOUT$IGRAPH_FUN == "with_fr") {
      # looks like default force-directed algorithm
      layout <- igraph::layout_with_fr(G)
    } else if (config$LAYOUT$IGRAPH_FUN == "sugiyama") {
      # minimzes edge crossings
      layout <- igraph::layout_with_sugiyama(G)$layout
    } else {
      stop("Type of igraph layout not supported")
    }
  }
  if (layout_type == "bipartite") {
    # layout <- igraph::layout_as_bipartite(
    #   graph = G,
    #   types = igraph::V(G)$vertex_types,
    #   hgap = config$LAYOUT$HGAP,
    #   vgap = config$LAYOUT$VGAP,
    #   maxiter = 100
    # )
    # layout <- layout[, 2:1] # horizontal to vertical
    n_nodes <- length(igraph::V(G))
    n_nodes2 <- n_nodes/2
    layout_emitter <- matrix(
      c(
        rep(0, times = n_nodes2),
        seq(
          from = 0,
          to = (n_nodes2-1)*config$LAYOUT$HGAP,
          by = config$LAYOUT$HGAP
        )
      ),
      nrow = n_nodes2,
      ncol = 2
    )
    layout_receiver <- matrix(
      c(
        rep(config$LAYOUT$VGAP, times = n_nodes2),
        seq(
          from = 0,
          to = (n_nodes2-1)*config$LAYOUT$HGAP,
          by = config$LAYOUT$HGAP
        )
      ),
      nrow = n_nodes2,
      ncol = 2
    )
    layout <- rbind(layout_emitter, layout_receiver)
    layout <- layout[igraph::vertex.attributes(G)$vertex_order, ]
    # if (config$LAYOUT$DISPERSE) {
    #   vgap <- config$LAYOUT$HGAP
    #   scale_factor <- 3
    #   new_layout[1:midpoint, 2][emitter_keys == 0] <- (
    #     new_layout[1:midpoint, 2][emitter_keys == 0] + scale_factor * vgap
    #   )
    #   new_layout[1:midpoint, 2][emitter_keys > 0] <- (
    #     new_layout[1:midpoint, 2][emitter_keys > 0] + 2 * scale_factor * vgap
    #   )
    #   new_layout[(midpoint + 1):num_v, 2][receiver_keys == 0] <- (
    #     new_layout[(midpoint + 1):num_v, 2][receiver_keys == 0]
    #     + scale_factor * vgap
    #   )
    #   new_layout[(midpoint + 1):num_v, 2][receiver_keys > 0] <- (
    #     new_layout[(midpoint + 1):num_v, 2][receiver_keys > 0]
    #     + 2 * scale_factor * vgap
    #   )
    # }
  }
  G$layout <- layout
  return(G)
}

build_visnetwork <- function(
  G,
  object_name,
  network_type,
  layout_type,
  config
) {
  nodes <- setDT(igraph::as_data_frame(G, what = "vertices"))
  edges <- setDT(igraph::as_data_frame(G, what = "edges"))
  network_components <- get_network_components(
    nodes = nodes,
    edges = edges,
    object_name = object_name,
    layout = G$layout,
    network_type = network_type,
    layout_type = layout_type,
    config = config,
    configure = FALSE
  )
  interactive_network <- do.call(apply_visnetwork, network_components)
}

get_network_components <- function(
  nodes,
  edges,
  object_name,
  layout,
  network_type,
  layout_type,
  config,
  configure
) {
  if (configure) {
    configure_component <- . %>% visNetwork::visConfigure(
      enabled = TRUE,
      showButton = TRUE
    )
  } else {
    configure_component <- NULL
  }
  network_components <- list(
    network_skeleton =  visNetwork::visNetwork(
      nodes = nodes,
      edges = edges,
      width = config$VISNETWORK$WIDTH,
      height = config$VISNETWORK$HEIGHT,
      main = edges$main_title[[1]],
      #submain = sprintf("%s", object_name),
      #footer = sprintf("Network type: %s", layout_type),
      background = config$VISNETWORK$BACKGROUND
    ),
    nodes_global = . %>% visNetwork::visNodes(
      shape = "dot",
      physics = FALSE,
      font = list(size = 18, align = "left")
    ),
    edges_global = . %>% visNetwork::visEdges(
      shadow = TRUE,
      arrows = "middle",
      smooth = list(enabled = TRUE, roundness = 0.75)
    ),
    layout = . %>% visNetwork::visIgraphLayout(
      layout = "layout.norm",
      layoutMatrix = layout
    ),
    options = get_visnetwork_options(selectionByName=TRUE),
    interactive = get_visnetwork_interactive(),
    configure = configure_component,
    legend = . %>% visNetwork::visLegend(
      enabled = TRUE,
      main = "Legend",
      position = "right",
      addEdges = edges_legend(
        network_type = network_type,
        config = config
      ),
      addNodes = nodes_legend(
        network_type = network_type,
        config = config
      ),
      useGroups = FALSE,
      zoom = FALSE
    ),
    physics = . %>% visNetwork::visPhysics(
      enabled = FALSE,
      maxVelocity = 10,
      timestep = 0.5,
      repulsion = list(damping=0.5)
    )
  )
  return(network_components)
}

apply_visnetwork <- function(
  network_skeleton,
  nodes_global = NULL,
  edges_global = NULL,
  layout = NULL,
  legend = NULL,
  options = NULL,
  interactive = NULL,
  configure = NULL,
  physics = NULL
) {
  vis_funcs <- list(
    nodes_global = nodes_global,
    edges_global = edges_global,
    layout = layout,
    legend = legend,
    options = options,
    interactive = interactive,
    configure = configure,
    physics = physics
  )
  # For NULL arguments, use the identity as pipeline step
  vis_funcs <- purrr::map_if(vis_funcs, is.null, ~ . %>% identity())
  return(
    network_skeleton %>%
      (vis_funcs$nodes_global) %>%
      (vis_funcs$edges_global) %>%
      (vis_funcs$layout) %>%
      (vis_funcs$legend) %>%
      (vis_funcs$options) %>%
      (vis_funcs$interactive) %>%
      (vis_funcs$configure) %>%
      (vis_funcs$physics)
  )
}

get_visnetwork_options <- function(
  selectionByName = FALSE,
  highlightNearestDegree = 1
) {
  if(selectionByName) {
    selectionOptions = list(
      variable = "name", multiple = TRUE,
      style = "width: 200px; height: 26px;
          background: #f8f8f8;
          color: darkblue;
          border:none;
          outline:none;"
    )
  } else {
    selectionOptions = NULL
  }
  return(
    . %>% visNetwork::visOptions(
      height = NULL, #"1000px",
      width = NULL, #"100%",
      highlightNearest = list(
        enabled = TRUE,
        #degree = highlightNearestDegree,
        degree = list(from = 1, to = 1),
        hover = TRUE,
        algorithm = "hierarchical"
      ),
      autoResize = TRUE,
      selectedBy = selectionOptions,
      collapse = TRUE,
      manipulation = TRUE
    )
  )
}

get_visnetwork_interactive <- function(
) {
  return(
    . %>% visNetwork::visInteraction(
      keyboard = list(enabled = TRUE),
      multiselect = TRUE,
      navigationButtons = FALSE
    )
  )
}

nodes_legend <- function(
  network_type,
  config
) {
  return(
    data.frame(
      label = c("Cell Type"),
      shape = c("dot"),
      color = c(config$NODE_COLORING$BACKGROUND)
    )
  )
}

edges_legend <- function(
  network_type,
  config
) {
  if (network_type == "ORA_network") {
    return(
      data.frame(
        color = c(
          config$EDGE_COLORING$ORA_COLOR_DOWN,
          config$EDGE_COLORING$ORA_COLOR_UP,
          config$EDGE_COLORING$ORA_COLOR_FLAT,
          config$EDGE_COLORING$ORA_COLOR_DIFF
        ),
        label = c("DOWN", "UP", "FLAT", "UP&DOWN"),
        arrows = c("to", "to", "to", "to")
      )
    )
  } else if (network_type %in% c(
    "condition1_network",
    "condition2_network"
    )) {
    return(
      data.frame(
        color = "darkblue",
        label = "#interactions",
        arrows = "to",
        width = 2
      )
    )
  } else if (network_type == "difference_network") {
    return(data.frame(
      color = "darkblue",
      label = "delta(interactions)",
      arrows = "to",
      width = 2,
      length = 10
    ))
  } else if (network_type %in% c(
    "up_regulated_network",
    "down_regulated_network"
  )) {
    return(data.frame(
      color = "darkblue",
      label = "toDo", # TODO
      arrows = "to",
      width = 2,
      length = 10
    ))
  } else {
    stop()
  }
}

edge_annotation_html <- function(
  edge_table,
  network_type
) {
  ORA_TYPE <- OR_UP <- BH_P_VALUE_UP <-
    OR_DOWN <- BH_P_VALUE_DOWN <-
    OR_FLAT <- BH_P_VALUE_FLAT <-
    NUM_LRIS_TOTAL <- NUM_LRIS_UP <-
    NUM_LRIS_DOWN <- NUM_LRIS_FLAT <- NULL
  if (network_type != "ORA_network") {
    edge_table[, "title" := ""]
  } else {
    edge_table[
      ,
      "title" := lapply(
        1:nrow(.SD),
        function(i) {
          h1 <- "<h4> ORA results: </h4>"
          if (ORA_TYPE[[i]] == "UP") {
            ora_results <- as.character(
              kableExtra::kbl(
                matrix(
                  c(
                    "Odds Ratio UP:", OR_UP[[i]],
                    "Adj. p-value UP:", BH_P_VALUE_UP[[i]]
                  ),
                  nrow = 2,
                  byrow = TRUE
                  )
              )
            )
          } else if (ORA_TYPE[[i]] == "DOWN") {
            ora_results <- as.character(
              kableExtra::kbl(
                matrix(
                  c(
                    "Odds Ratio DOWN:", OR_DOWN[[i]],
                    "Adj. p-value DOWN:", BH_P_VALUE_DOWN[[i]]
                  ),
                  nrow = 2,
                  byrow = TRUE
                )
              )
            )
          } else if (ORA_TYPE[[i]] == "FLAT") {
            ora_results <- as.character(
              kableExtra::kbl(
                matrix(
                  c(
                    "Odds Ratio FLAT:", OR_FLAT[[i]],
                    "Adj. p-value FLAT:", BH_P_VALUE_FLAT[[i]]
                  ),
                  nrow = 2,
                  byrow = TRUE
                )
              )
            )
          } else {
            ora_results <- as.character(
              kableExtra::kbl(
                matrix(
                  c(
                    "Odds Ratio UP:", OR_UP[[i]],
                    "Adj. p-value UP:", BH_P_VALUE_UP[[i]],
                    "Odds Ratio DOWN:", OR_DOWN[[i]],
                    "Adj. p-value DOWN:", BH_P_VALUE_DOWN[[i]]
                  ),
                  nrow = 4,
                  byrow = TRUE
                )
              )
            )
          }
          n_inter <- as.character(
            kableExtra::kbl(
              matrix(
                c(
                  "TOTAL:", NUM_LRIS_TOTAL[[i]],
                  "UP:", NUM_LRIS_UP[[i]],
                  "DOWN:", NUM_LRIS_DOWN[[i]],
                  "FLAT:", NUM_LRIS_FLAT[[i]]
                ),
                nrow = 4,
                byrow = TRUE
              )
            )
          )
          paste0(
            #h1,
            ora_results,
            "<br> Number of interactions: <br>",
            n_inter
          )
        }
      )
    ]
  }
  return(edge_table)
}

# get_cci_change_graph <- function(cci_table_detected) {
#   LOGFC <- value <- label <- color <- REGULATION <-
#     EMITTER_CELLTYPE <- RECEIVER_CELLTYPE <- LRI <- NULL
#   dt_edge = cci_table_detected[REGULATION %in% c('UP', 'DOWN'),
#                          list(
#                            EMITTER_CELLTYPE, RECEIVER_CELLTYPE,
#                            LRI, LOGFC,
#                            REGULATION
#                          )
#   ][,
#     value := abs(LOGFC)][
#       ,
#       label := format(round(LOGFC, 2), nsmall=2)][
#         ,
#         color := ifelse(REGULATION == 'UP', 'red', 'blue')
#       ]
#   cci_change_graph = igraph::graph_from_data_frame(dt_edge, directed = TRUE, vertices = NULL)
#   return(cci_change_graph)
# }
#
# get_cci_change_subgraph_list <- function(cci_table_detected, LRIs = NULL) {
#   LRI <- NULL
#   if(is.null(LRIs)) {
#     lris = unique(cci_table_detected[, LRI])
#   } else {
#     lris = LRIs
#   }
#   subgraphs = purrr::map(
#     lris,
#     ~ get_cci_change_subgraph(cci_table_detected, .x)
#   )
#   if(length(lris)==1) {
#     return(subgraphs[[1]])
#   } else {
#     return(subgraphs)
#   }
# }
#
# LRI_subgraph_metrics <- function(cci_table_detected, LRIs=NULL) {
#   subgraphs = get_cci_change_subgraph_list(cci_table_detected, LRIs)
#
#   subg_metrics_l = purrr::map(
#     subgraphs,
#     function(subg) {
#
#       metrics = list(
#         LRI = subg$LRI,
#         num_edges = length(igraph::E(subg)),
#         num_loops = sum(igraph::which_loop(subg)),
#         num_reciprocal = sum(igraph::which_mutual(subg)) - sum(igraph::which_loop(subg)),
#         motifs_3_total = igraph::count_motifs(subg, size=3)
#       )
#       motifs_3 = as.list(igraph::motifs(subg, size = 3))
#       names(motifs_3) = paste0('Motifs_3_id_', 0:15)
#
#       return(c(metrics, motifs_3))
#
#     }
#   )
#
#   subg_metrics_dt = data.table::rbindlist(subg_metrics_l)
#   return(subg_metrics_dt)
# }
#
# get_cci_change_subgraph <- function(
#   cci_table_detected,
#   LRIs
# ) {
#   SUBG_LIMIT <- 5
#   if(length(LRIs) > SUBG_LIMIT) {
#     stop(paste0("The number of LRIs must be less than ", SUBG_LIMIT))
#   }
#   dt_edge = cci_table_detected[
#     (REGULATION %in% c('UP', 'DOWN'))
#     & (LRI %in% LRIs),
#     c(
#       "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",
#       "LRI", "LOGFC_ABS",
#       "REGULATION"
#     )
#   ][,
#     value := LOGFC_ABS][
#       ,
#       # label := format(round(LOGFC, 2), nsmall=2)][
#       label := LRI][
#         ,
#         color := ifelse(REGULATION == 'UP', 'red', 'blue')
#       ]
#   if(nrow(dt_edge) == 0) {
#     return(igraph::make_empty_graph())
#   } else {
#     g = igraph::graph_from_data_frame(dt_edge, directed = TRUE, vertices = NULL)
#     g$LRI = LRIs
#     return(g)
#   }
# }






