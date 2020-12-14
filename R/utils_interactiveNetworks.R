# 
construct_graph <- function(object) {
  
  na_handler <- function(dt, message) {
    has_na = sum(is.na(dt)) > 0
    if (has_na) {
      message(paste0(
        message,
        "na_handler: #NAs", ":", sum(is.na(dt)), "; NA->1"
      ))
      dt[is.na(dt)] <- 1
    }
    return(dt)
  }
  inf_handler <- function(dt, message) {
    has_inf = sum(dt == Inf) > 0
    if (has_inf) {
      stop("construct_graph: Inf values in dt_edge")
    }
    return(dt)
  }
  classify_edges <- function(dt_edge,
                             config = setup_graph_config(),
                             use_adjpval = TRUE) {
    # TODO: Use setup_graph_config()
    LABEL_DIFF <- "diff"
    LABEL_UP <- "up"
    LABEL_DOWN <- "down"
    LABEL_ROBUST <- "robust"
    LABEL_NONE <- "none"
    
    # TODO
    CUTOFF_CHANGE_DIFF <- config$EDGE_COLORING$CUTOFF_DIFF
    CUTOFF_CHANGE_UP <- config$EDGE_COLORING$CUTOFF_UP  # Change name
    CUTOFF_CHANGE_DOWN <- config$EDGE_COLORING$CUTOFF_DOWN  # Change name
    CUTOFF_ROBUST_UP <- 0.99
    CUTOFF_ROBUST_DOWN <- 0.99
    
    edges_classification <- purrr::pmap_chr(
      dt_edge,
      function(OR_DIFF, OR_UP, OR_DOWN,
               P_VALUE_DIFF, P_VALUE_UP, P_VALUE_DOWN,
               BH_P_VALUE_DIFF, BH_P_VALUE_UP, BH_P_VALUE_DOWN,
               ...) {
        if (use_adjpval) {
          is_P_VALUE_DIFF <- BH_P_VALUE_DIFF < 0.05
          is_P_VALUE_UP <- BH_P_VALUE_UP < 0.05
          is_P_VALUE_DOWN <- BH_P_VALUE_DOWN < 0.05
        } else {
          is_P_VALUE_DIFF <- P_VALUE_DIFF < 0.05
          is_P_VALUE_UP <- P_VALUE_UP < 0.05
          is_P_VALUE_DOWN <- P_VALUE_DOWN < 0.05
        }
        
        is_OR_DIFF_CHANGE <- OR_DIFF >= CUTOFF_CHANGE_DIFF
        is_OR_UP_CHANGE <- OR_UP >= CUTOFF_CHANGE_UP
        is_OR_DOWN_CHANGE <- OR_DOWN >= CUTOFF_CHANGE_DOWN
        is_OR_UP_ROBUST <- OR_UP <= CUTOFF_ROBUST_UP
        is_OR_DOWN_ROBUST <- OR_DOWN <= CUTOFF_ROBUST_DOWN
        
        if (is_OR_UP_CHANGE & is_P_VALUE_UP & is_OR_DOWN_CHANGE & is_P_VALUE_DOWN) {
          return(LABEL_DIFF)
        } else if (is_OR_UP_CHANGE & is_P_VALUE_UP) {
          return(LABEL_UP)
        } else if (is_OR_DOWN_CHANGE & is_P_VALUE_DOWN) {
          return(LABEL_DOWN)
        } else if (is_OR_UP_ROBUST & is_P_VALUE_UP & is_OR_DOWN_ROBUST & is_P_VALUE_DOWN) {
          return(LABEL_ROBUST)
        } else {
          return(LABEL_NONE)
        }
      }
    )
  }
  
  dt_edge <- process_celltype_pairs_enrichment(object@ora_default$ER_CELLTYPES)
  
  dt_edge[,
    "Ligand_cell_clean" := Ligand_cell][,
    "Receptor_cell_clean" := Receptor_cell][,
    "Ligand_cell" := paste0(Ligand_cell, " (L)")][,
    "Receptor_cell" := paste0(Receptor_cell, " (R)")
  ]
  dt_edge = na_handler(dt_edge, 'construct_graph: ')
  dt_edge = inf_handler(dt_edge, 'construct_graph: ')

  dt_edge[, edge_type := classify_edges(dt_edge)]
  
  add_celltype_average_counts <- function(dt_edge, object) {
    counts = average_celltype_counts(object)
    dt_edge = merge(dt_edge, counts, 
                    by.x = 'Ligand_cell_clean',
                    by.y = 'name',
                    all = FALSE)
    dt_edge[, 'Ligand_cell_counts' := counts]
    dt_edge = merge(dt_edge, counts,
                    by.x = 'Receptor_cell_clean',
                    by.y = 'name',
                    all = FALSE)
    dt_edge[, 'Receptor_cell_counts' := counts.y]
    dt_edge[, counts.x := NULL]
    dt_edge[, counts.y := NULL]
    setcolorder(dt_edge, c(3:dim(dt_edge)[2], 2, 1))
    return(dt_edge)
  }
  dt_edge = add_celltype_average_counts(dt_edge, object)
  
  add_celltype_pair_counts <- function(dt_edge, object) {
    num_interactions_cellpairs <- function(object) {
      # cond1 = object@parameters$cond1_name
      # cond2 = object@parameters$cond2_name
      # cond1 = 'YOUNG'
      # colname1 = sprintf("CCI_DETECTED_%s", cond1)
      # cond2 = 'OLD'
      # colname2 = sprintf("CCI_DETECTED_%s", cond2)
      
      counts_y = object@cci_detected[CCI_DETECTED_YOUNG == TRUE,
                                     .(count = .N),
                                     by = .(EMITTER_CELLTYPE, RECEIVER_CELLTYPE)
      ][, .(
        "Ligand_cell_clean" = EMITTER_CELLTYPE,  # TODO: Change these colnames
        "Receptor_cell_clean" = RECEIVER_CELLTYPE,
        "num_interacts_young" = count
      )]
      counts_o = object@cci_detected[CCI_DETECTED_OLD == TRUE,
                                     .(count = .N),
                                     by = .(EMITTER_CELLTYPE, RECEIVER_CELLTYPE)
      ][, .(
        "Ligand_cell_clean" = EMITTER_CELLTYPE,  # TODO: Change these colnames
        "Receptor_cell_clean" = RECEIVER_CELLTYPE,
        "num_interacts_old" = count
      )]
      counts = merge(counts_y, counts_o, 
                     by = c('Ligand_cell_clean', 'Receptor_cell_clean'),
                     all = TRUE)
      return(counts)
    }
    interacts = num_interactions_cellpairs(object)
    if (nrow(dt_edge) != nrow(interacts)) {
      stop('Must have same nr rows.')  # TODO
    }
    dt_merge = merge(dt_edge, interacts, 
                     by=c('Ligand_cell_clean', 'Receptor_cell_clean'),
                     all=FALSE)  # TODO: Check that same interactions are covered in both dts
    first_2_cols = c('Ligand_cell', 'Receptor_cell')
    data.table::setcolorder(dt_merge, c(first_2_cols, setdiff(names(dt_merge), first_2_cols)))
    return(dt_merge)
  }
  dt_edge = add_celltype_pair_counts(dt_edge, object)
  dt_edge[, num_interacts_diff := num_interacts_old - num_interacts_young]
  
  # TODO: Enrich with other info: cell families, GOs, ORA on cell types.
  dt_vertex = funion(
    dt_edge[, .(
      name=Ligand_cell,
      counts=Ligand_cell_counts)],
    dt_edge[, .(
      name=Receptor_cell,
      counts=Receptor_cell_counts)]
  )
  
  
  G <- igraph::graph_from_data_frame(
    d = dt_edge,
    directed = TRUE,
    vertices = dt_vertex
  )

  G$scDiffCom <- object
  
  # TODO: Remove, this info is redundant.
  # extracted <- extract_from_object(object)
  # tissue = extracted$tissue
  # G$interacts_counts <- num_interactions_cellpairs(object)
  # G$celltype_counts <- average_celltype_counts(object)
  return(G)
}
process_celltype_pairs_enrichment <- function(ora_ER_cells) {
  
  COL_LIGAND_RECEPTOR_CELLTYPES <- "ER_CELLTYPES"
  COL_P_VALUE_DIFF <- "P_VALUE_DIFF"
  COL_P_VALUE_FLAT <- "P_VALUE_FLAT"
  COL_P_VALUE_UP <- "P_VALUE_UP"
  COL_P_VALUE_DOWN <- "P_VALUE_DOWN"
  
  COL_BH_P_VALUE_DIFF <- "BH_P_VALUE_DIFF"
  COL_BH_P_VALUE_FLAT <- "BH_P_VALUE_FLAT"
  COL_BH_P_VALUE_UP <- "BH_P_VALUE_UP"
  COL_BH_P_VALUE_DOWN <- "BH_P_VALUE_DOWN"
  
  dt_ctypes = copy(ora_ER_cells)
  dt_ctypes <- dt_ctypes[, c("Ligand_cell", "Receptor_cell") := list(
    sub("_.*", "", VALUE),
    sub(".*_", "", VALUE)
  )]
  
  dt_ctypes <- dt_ctypes[,
    'COUNTS_UP' := COUNTS_VALUE_REGULATED_UP][,
    'COUNTS_DOWN' := COUNTS_VALUE_REGULATED_DOWN][,
    'COUNTS_DIFF' := COUNTS_VALUE_REGULATED_DIFF][,
    'COUNTS_FLAT' := COUNTS_VALUE_REGULATED_FLAT][,
    'COUNTS_TOTAL' := COUNTS_VALUE_REGULATED_DIFF + COUNTS_VALUE_REGULATED_FLAT  # Size of union{LR young and old}
  ]
  
  subset_celltypes_dt_columns <- function(dt, cols) {
    if (length(cols) != 19) {
      stop('Not implemented')
    }
    dt <- dt[, list(
      get(cols[1]), get(cols[2]), get(cols[3]),
      get(cols[4]), get(cols[5]), get(cols[6]),
      get(cols[7]), get(cols[8]), get(cols[9]),
      get(cols[10]), get(cols[11]), get(cols[12]),
      get(cols[13]), get(cols[14]), get(cols[15]),
      get(cols[16]), get(cols[17]), get(cols[18]),
      get(cols[19])
    )]
    names(dt) = cols
    return(dt)
  }
  
  cols_to_select <- c(
    "Ligand_cell", "Receptor_cell", "OR_DIFF", "OR_UP", "OR_DOWN", "OR_FLAT",
    COL_P_VALUE_DIFF, COL_P_VALUE_UP, COL_P_VALUE_DOWN, COL_P_VALUE_FLAT,
    COL_BH_P_VALUE_DIFF, COL_BH_P_VALUE_UP, COL_BH_P_VALUE_DOWN, COL_BH_P_VALUE_FLAT,
    'COUNTS_UP', 'COUNTS_DOWN', 'COUNTS_DIFF', 'COUNTS_FLAT', 'COUNTS_TOTAL'
  )
  dt_ctypes = subset_celltypes_dt_columns(dt_ctypes, cols_to_select)
  return(dt_ctypes)
}
map_edgetype_to_color <- function(edge_label, config = setup_graph_config()) {
  
  # TODO: Inline variables
  COLOR_DIFF <- config$EDGE_COLORING$COLOR_DIFF
  COLOR_UP <- config$EDGE_COLORING$COLOR_UP
  COLOR_DOWN <- config$EDGE_COLORING$COLOR_DOWN
  COLOR_ROBUST <- config$EDGE_COLORING$COLOR_ROBUST
  COLOR_NONE <- config$EDGE_COLORING$COLOR_NONE
  
  label_to_color_map <- list(
    "diff" = COLOR_DIFF,
    "up" = COLOR_UP,
    "down" = COLOR_DOWN,
    "robust" = COLOR_ROBUST,
    "none" = COLOR_NONE
  )
  return(label_to_color_map[[edge_label]])
}
setup_graph <- function(G, config = setup_graph_config(), use_adjpval, disperse) {
  
  setup_vertices <- function(G, config) {
    
    infer_vertex_types <- function(vertex_names) {
      vertex_types <- purrr::map_lgl(
        strsplit(vertex_names, "[()]"),
        ~ .x[2] == "L"
      )
      return(vertex_types)
    }
    add_vertex_size <- function(G) {
      MAXSIZE <- 20
      MINSIZE <- 5
      igraph::V(G)$vertex.size <- rescale(igraph::V(G)$counts, MINSIZE, MAXSIZE)
      return(G)
    }
    
    igraph::V(G)$vertex_types <- infer_vertex_types(igraph::V(G)$name)
    G <- add_vertex_size(G)
    return(G)
  }
  setup_edges <- function(G, config, use_adjpval) {
    
    add_edge_width <- function(G) {
      WIDTH_MIN <- 0
      WIDTH_MAX <- 10
      num_interacts <- igraph::as_data_frame(G, 'edges')$COUNTS_TOTAL
      igraph::E(G)$width <- rescale(sqrt(num_interacts), WIDTH_MIN, WIDTH_MAX)
      return(G)
    }
    
    igraph::E(G)$arrow.size <- config$EDGE_STYLE$ARROW_SIZE
    igraph::E(G)$edge.color <- purrr::map_chr(igraph::E(G)$edge_type, map_edgetype_to_color)
    G <- add_edge_width(G)
    return(G)
  }
  setup_layout <- function(G, config, disperse) {
    layout <- igraph::layout_as_bipartite(
      G,
      types = igraph::V(G)$vertex_types,
      hgap = config$LAYOUT$HGAP,
      vgap = config$LAYOUT$VGAP
    )
    layout <- layout[, 2:1] # horizontal to vertical
    layout <- sort_layout_vertices(layout, G, config, disperse)
    G$layout <- layout
    return(G)
  }
  
  G <- setup_vertices(G, config)
  G <- setup_edges(G, config, use_adjpval)
  G <- setup_layout(G, config, disperse)
  return(G)
}

average_celltype_counts <- function(object) {
  
  dt_t <- data.table::copy(object@cci_detected)

  counts_dt <- data.table::funion(
    dt_t[, .(
      name = EMITTER_CELLTYPE,
      counts = (EMITTER_NCELLS_OLD + EMITTER_NCELLS_YOUNG) / 2 # average counts in young-old
    )],
    dt_t[, .(
      name = RECEIVER_CELLTYPE,
      counts = (RECEIVER_NCELLS_OLD + RECEIVER_NCELLS_YOUNG) / 2 # average counts in young-old
    )],
  )
  counts_dt <- counts_dt[order(counts_dt$name)]
  return(counts_dt)
}
num_interactions_object <- function(object) {
  TOTAL = object@cci_detected[, .N]
  UP = object@cci_detected[REGULATION_SIMPLE=='UP', .N]
  DOWN = object@cci_detected[REGULATION_SIMPLE=='DOWN', .N]
  DIFF = object@cci_detected[REGULATION_SIMPLE=='DOWN' | REGULATION_SIMPLE=='UP', .N]
  FLAT = object@cci_detected[REGULATION_SIMPLE=='FLAT', .N]
  return(list(
    TOTAL = TOTAL,
    UP = UP,
    DOWN = DOWN,
    DIFF = DIFF,
    FLAT = FLAT
  ))
}
rescale <- function(v, min_, max_) {
  (v - min(v)) / (max(v) - min(v)) * (max_ - min_) + min_
}

#### Functions for improved layout of biparite type ####
sort_layout_vertices <- function(layout, G, config, disperse) {

  get_edges_from_vertex <- function(v_name, G) {
    return(igraph::incident(G, v_name, mode = "out"))
  }
  get_edges_to_vertex <- function(v_name, G) {
    return(igraph::incident(G, v_name, mode = "in"))
  }
  count_up_and_down_edges <- function(edges, config) {
    color_up <- config$EDGE_COLORING$COLOR_UP
    color_down <- config$EDGE_COLORING$COLOR_DOWN
    
    counts <- table(edges$edge.color)
    num_up <- as.integer(counts[color_up])
    num_down <- as.integer(counts[color_down])
    
    num_up <- ifelse(is.na(num_up), 0, num_up)
    num_down <- ifelse(is.na(num_down), 0, num_down)
    
    return(list(N_UP = num_up, N_DOWN = num_down))
  }
  vertex_sort_key_atomic <- function(v_name, from, G, config) {
    if (from) {
      edges <- get_edges_from_vertex(v_name, G)
    }
    else {
      edges <- get_edges_to_vertex(v_name, G)
    }
    counts <- count_up_and_down_edges(edges, config)
    sort_key <- as.integer(counts$N_UP - counts$N_DOWN)
    return(sort_key)
  }
  vertex_sort_keys <- function(v_names, from, G, config) {
    sort_keys <- purrr::map_dbl(
      v_names,
      ~ vertex_sort_key_atomic(.x, from, G, config)
    )
    return(sort_keys)
  }
  
  vertex_names <- igraph::V(G)$name
  num_v <- length(vertex_names)
  midpoint <- num_v / 2
  ligand_vertices <- vertex_names[1:midpoint]
  receptor_vertices <- vertex_names[(midpoint + 1):num_v]

  ligand_keys <- vertex_sort_keys(ligand_vertices, from = TRUE, G, config)
  receptor_keys <- vertex_sort_keys(receptor_vertices, from = FALSE, G, config)

  # Reorders vertices to sort
  layout_copy <- copy(layout)
  new_layout <- copy(layout)
  new_layout[order(ligand_keys), ] <- layout_copy[1:midpoint, ]
  new_layout[midpoint + order(receptor_keys), ] <- layout_copy[(midpoint + 1):num_v, ]

  if (disperse) {

    vgap <- config$LAYOUT$HGAP
    scale_factor <- 3

    new_layout[1:midpoint, 2][ligand_keys == 0] <- (
      new_layout[1:midpoint, 2][ligand_keys == 0] + scale_factor * vgap
    )
    new_layout[1:midpoint, 2][ligand_keys > 0] <- (
      new_layout[1:midpoint, 2][ligand_keys > 0] + 2 * scale_factor * vgap
    )

    new_layout[(midpoint + 1):num_v, 2][receptor_keys == 0] <- (
      new_layout[(midpoint + 1):num_v, 2][receptor_keys == 0]
      + scale_factor * vgap
    )
    new_layout[(midpoint + 1):num_v, 2][receptor_keys > 0] <- (
      new_layout[(midpoint + 1):num_v, 2][receptor_keys > 0]
      + 2 * scale_factor * vgap
    )
  }
  return(new_layout)
}

#### Interactive networks ####
extract_from_object <- function(object, class_signature, subobject_name) {
  tryCatch(
    error = function(cnd) {
      if (class_signature == "scDiffCom") {
        return(list(
          tissue = object@parameters$object_name,
          ora = object@ora_default,
          dt_filtered = object@cci_detected
        ))
      }
      if (class_signature == "scDiffComCombined") {
        ora <- object@ora_default
        ora$ER_CELLTYPES <- ora$ER_CELLTYPES[ID == subobject_name][, ID := NULL]
        return(list(
          tissue = subobject_name,
          ora = ora,
          dt_filtered = object@cci_detected[ID == subobject_name][, ID := NULL]
        ))
      }
    },
    stop("Can't extract necessary data from object. Maybe check scDiffCom object.")
  )
}
check_network_layout_type_arg <- function(network_layout_type) {
  if (!(network_layout_type %in% c("bipartite", "celltypes"))) {
    message <- sprintf('network_layout_type "%s" not in "bipartite", "celltypes"', network_layout_type)
    stop(message)
  }
}
remove_LR_label <- function(s) {
  return(
    sub(" [(][LR][)]", "", s)
  )
}
map_bipartite_to_community <- function(G) {
  nodes <- data.table::data.table(igraph::as_data_frame(G, what = "vertices"))
  edges <- data.table::data.table(igraph::as_data_frame(G, what = "edges"))

  # Process edges
  edges[, 
        from := remove_LR_label(from)][,
        to := remove_LR_label(to)
  ]

  # Process nodes
  nodes[, name := remove_LR_label(name)]
  nodes[, vertex_types := NULL]
  nodes <- nodes[, list(
    vertex.size = mean(vertex.size),
    counts = mean(counts)  # TODO: Not sure still works
  ), by = name]

  G_new <- igraph::graph_from_data_frame(vertices = nodes, d = edges)
  # G_new$tissue <- G$tissue
  # G_new$ora <- G$ora
  # G_new$dt <- G$dt
  # G_new$interacts_counts <- G$interacts_counts
  # G_new$celltype_counts <- G$celltype_counts

  # New layout
  # layout = igraph::layout_nicely(G_new, dim=2)  # determines best layout, likely calls fr
  # layout = igraph::layout_with_fr(G_new)  # looks like default force-directed algorithm
  layout <- igraph::layout_with_sugiyama(G_new)$layout # minimzes edge crossings
  # layout = igraph::layout_in_circle(G_new)
  G_new$layout <- layout

  return(G_new)
}
build_igraph <- function(object, network_layout_type) {
  G <- construct_graph(object)
  G <- setup_graph(G, use_adjpval = TRUE, disperse = TRUE)
  if (network_layout_type == "celltypes") {
    G <- map_bipartite_to_community(G)
  }
  return(G)
}
extract_node_edges <- function(G, exclude_nonsign) {
  nodes <- data.table::data.table(igraph::as_data_frame(G, what = "vertices"))
  edges <- data.table::data.table(igraph::as_data_frame(G, what = "edges"))

  if (exclude_nonsign) {
    # TODO
    edges <- edges[edge_type != "none"]
  }

  return(
    list(nodes = nodes, edges = edges)
  )
}
build_network_skeleton <- function(object,
                                   nodes,  # replace with igraph? would carry more attributes
                                   edges,
                                   network_representation_type,
                                   network_layout_type,
                                   config = setup_graph_config()) {
  nodes_ <- data.table::copy(nodes)
  edges_ <- data.table::copy(edges)
  
  # TODO: Extract function
  node_annotation_html <- function(nodes) {
    annotations = purrr::pmap_chr(
      nodes,
      function(name, counts, ...) {
        temp_dt = data.table::data.table(
          Type=c("Cell counts", "Cell OR"),
          Total=NA,
          UP=NA,
          DOWN=NA,
          FLAT=NA,
          DIFF=NA
        )
        header = sprintf("<h3> %s </h3>", name)
        table_html = as.character(kableExtra::kbl(temp_dt))
        num_cells_paragraph = sprintf("<p> Average num cells: %g", counts)
        html = paste0(header, num_cells_paragraph)  # table_html
        return(html) 
      }
    )
    return(annotations)
  }
  nodes_[
    ,
    "id" := name
  ][
    ,
    "label" := name
  ][
    ,
    "value" := vertex.size
  ][
    ,
    c("color.background", "color.border", 
      "color.highlight.background", "color.highlight.border",
      "color.hover.background", "color.hover.border") := list(
        config$NODE_COLORING$BACKGROUND, config$NODE_COLORING$BORDER,
        config$NODE_COLORING$HIGHLIGHT$BACKGROUND, config$NODE_COLORING$HIGHLIGHT$BORDER,
        config$NODE_COLORING$HOVER$BACKGROUND, config$NODE_COLORING$HOVER$BORDER)
  ][
    ,
    "shadow" := TRUE
  ][
    ,
    # 'shapeProperties' := NULL][,
    "title" := node_annotation_html(nodes_)
  ]

  if (network_layout_type == "bipartite") {
    nodes_[
      ,
      "group" := vertex_types
    ][
      ,
      "level" := ifelse(vertex_types == TRUE, 1, 2) # For hierarchical layout from vis
    ]
  } else if (network_layout_type == "celltypes") {
  } else {
    stop('argument network_layout_type not in "bipartite", "celltypes"')
  }
  
  num_interacts = num_interactions_object(object)
  edge_annotation_html <- function(edges, network_representation_type) {
    edges_tmp = copy(edges)
    edges_tmp[, network_representation_type := network_representation_type]
    annotations = purrr::pmap(
      edges_tmp,
      function(COUNTS_TOTAL, COUNTS_UP, COUNTS_DOWN, COUNTS_DIFF, COUNTS_FLAT,
               OR_UP, OR_DOWN, OR_DIFF, OR_FLAT,
               BH_P_VALUE_UP, BH_P_VALUE_DOWN, BH_P_VALUE_DIFF, BH_P_VALUE_FLAT,
               Ligand_cell_clean, Receptor_cell_clean, 
               TISSUE_ENRICHED_LR, 
               network_representation_type,
               ...) {
        
        if ( network_representation_type != "ORA" ) {
          return(NULL)
        }
        
        format_float = function(x) sprintf("%.1f", x)
        format_pval = function(x) round(x, digits=3)
        cellpair_counts_ora = data.table::data.table(
          TYPE = c("Tissue counts", "Counts", "OR", "BH_P_VALUE"),
          TOTAL = c(num_interacts$TOTAL, as.integer(COUNTS_TOTAL), NA, NA),
          UP = c(num_interacts$UP, as.integer(COUNTS_UP), format_float(OR_UP), format_pval(BH_P_VALUE_UP)),
          DOWN = c(num_interacts$DOWN, as.integer(COUNTS_DOWN), format_float(OR_DOWN), format_pval(BH_P_VALUE_DOWN)),
          DIFF = c(num_interacts$DIFF, as.integer(COUNTS_DIFF), format_float(OR_DIFF), format_pval(BH_P_VALUE_DIFF)),
          FLAT = c(num_interacts$FLAT, as.integer(COUNTS_FLAT), format_float(OR_FLAT), format_pval(BH_P_VALUE_FLAT))
        )
        
        header = sprintf("<h3> %s -(to)-> %s </h3>", Ligand_cell_clean, Receptor_cell_clean)
        
        # TOP LR for cellpair
        # TODO: Extract function. Adds the LR table
        NUM_TOP_LR = 5
        ER_celltype = paste0(Ligand_cell_clean, "_", Receptor_cell_clean)
        lr_genes_significant_sorted = object@cci_detected[ER_CELLTYPES == ER_celltype 
                                          & REGULATION_SIMPLE %in% c('UP', 'DOWN')][
                                            order(-LOGFC_ABS),
                                            .(LR_GENES, LOGFC, REGULATION)][  # ORA_SCORE might overflow
                                              1:min(NUM_TOP_LR, data.table::.N),]  
        top_LR_table = kableExtra::kbl(lr_genes_significant_sorted, caption=sprintf('Top %s LR by LOGFC_ABS', NUM_TOP_LR))
        
        # Activity of LR found in most cellpairs of tissue
        # TODO: Extract function
        NUM_FREQ_LR = 5
        LR_freq = object@ora_default$LR_GENES[order(-OR_DIFF), VALUE][1:NUM_FREQ_LR]  # TODO: Counts or OR?
        dt = object@cci_detected[
          ER_CELLTYPES == ER_celltype][
            match(LR_freq, LR_GENES)][
              , list(LR_GENES, LOGFC, REGULATION)]
        dt[, LR_GENES := LR_freq]
        # TODO: -> log base 2
        dt[, LOGFC := as.character(LOGFC)]
        dt[is.na(LOGFC), LOGFC := 'Not detected']
        freq_LR_table = as.character(kableExtra::kbl(dt, caption=sprintf('Most frequent %s LR that change in tissue', NUM_FREQ_LR)))
        
        # Activity of GO soecufuc to the cell pair
        # TODO: Extract function
        NUM_TOP_GO = 5
        ora_dt = build_ora_dt(object@cci_detected[ER_CELLTYPES == ER_celltype], 
                              object@parameters$threshold_logfc, 
                              c('UP', 'DOWN', 'DIFF', 'FLAT'),
                              'GO_TERMS',
                              'mouse')
        ora_dt = ora_dt[, c('DIFFERENTIAL') := ifelse(
          BH_P_VALUE_UP < 0.05 & BH_P_VALUE_DOWN < 0.05, 'DIFF', ifelse(
            BH_P_VALUE_UP < 0.05, 'UP', ifelse(
              BH_P_VALUE_DOWN < 0.05, 'DOWN', 'FLAT')
            )
          )
        ]
        ora_dt = ora_dt[DIFFERENTIAL %in% c('UP', 'DOWN', 'DIFF')]
        ora_dt = ora_dt[order(-OR_DIFF)]
        num_rows_display = min(dim(ora_dt)[1], NUM_TOP_GO)
        ora_dt = ora_dt[1:num_rows_display]
        
        GO_table = ora_dt[, list(VALUE, DIFFERENTIAL)]  # , OR_UP, BH_P_VALUE_UP, OR_DOWN, BH_P_VALUE_DOWN
        GO_table = kableExtra::kbl(GO_table, caption=sprintf('TOP %s enriched GO for ER tuple', NUM_TOP_GO))
        
        # TODO: Extract function
        html_break = '</br>'
        html = paste0(header,
                      as.character(kableExtra::kbl(cellpair_counts_ora, caption='Cell pair counts and ORA')),
                      html_break,
                      top_LR_table,
                      html_break,
                      freq_LR_table,
                      html_break,
                      GO_table)
        return(html)
      }
    )
    return(annotations)
  }
  edge_color <- function(edges, network_representation_type) {
    if (network_representation_type == 'ORA') {
      return(edges[, edge.color])
    } else {
      map2colors <- function(x,n){
        RColorBrewer::brewer.pal(n, name='RdYlGn')[cut(x,n)]
      }
      # Discretizing and mapping to colors
      values = edge_value(edges, network_representation_type)
      colors = map2colors(values, n=3)
      return(colors)
    }
  }
  edge_value <- function(edges, network_representation_type) {
    if (network_representation_type == 'ORA') {
      return(NULL)
    } else if (network_representation_type == 'COUNTS_YOUNG') {
      return(edges[, num_interacts_young])
    } else if (network_representation_type == 'COUNTS_OLD') {
      return(edges[, num_interacts_young])
    } else if (network_representation_type == 'COUNTS_DIFF') {
      num_diff = edges[, num_interacts_diff]
      num_diff[num_diff == 0] = NA
      return(num_diff)
    } else {
      stop()
    }
  }
  edge_label <- function(edges, network_representation_type) {
    if (network_representation_type == 'ORA') {
      return(NULL)
    } else {
      value = edge_value(edges, network_representation_type)
      return(as.character(value))
    }
  }
  abs_or_null <- function(v) {
    if (is.null(v)) {return(NULL)}
    v_ = copy(v)
    v_[!is.null(v_)] = abs(v_[!is.null(v_)])
    return(v_)
  }
  edges_[
    ,
    value := abs_or_null(edge_value(edges_, network_representation_type))  # TODO: bad practice
  ][
    ,
    label := as.character(edge_value(edges_, network_representation_type))
  ][
    ,
    color.color := edge_color(edges_, network_representation_type)
  ][
    ,
    color.highlight := edge_color(edges_, network_representation_type)
  ][
    ,
    color.hover := edge_color(edges_, network_representation_type)
  ][
    ,
    smooth := TRUE
  ][
    ,
    title := edge_annotation_html(edges_, network_representation_type)
  ]
  print(edges_)
  return(
    visNetwork::visNetwork(
      nodes_, edges_,
      width = "100%", # height='100%',
      main = sprintf("Network representation of: %s", network_representation_type),
      submain = sprintf("%s", object@parameters$object_name),
      footer = sprintf("Network type: %s", network_layout_type),
      background = "white"
    )
  )
}
get_network_components <- function(object,
                                   nodes,
                                   edges,
                                   layout,
                                   network_representation_type,
                                   network_layout_type,
                                   configure = FALSE) {
  if (configure) {
    configure_component <- . %>% visNetwork::visConfigure(
      enabled = TRUE,
      showButton = TRUE
    )
  } else {
    configure_component <- NULL
  }

  network_components <- list(
    network_skeleton = build_network_skeleton(
      object,
      nodes, edges,
      network_representation_type = network_representation_type,
      network_layout_type = network_layout_type
    ),
    nodes_global = . %>% visNetwork::visNodes(
      shape = "dot",
      physics = TRUE  # TODO
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
    options = . %>% visNetwork::visOptions(
      height = '1000px',
      width = '100%',
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      selectedBy = list(
        variable = "name", multiple = TRUE, 
        style = 'width: 200px; height: 26px; 
        background: #f8f8f8;
        color: darkblue;
        border:none;
        outline:none;'),
      collapse = TRUE,
      manipulation = TRUE
    ),
    interactive = . %>% visNetwork::visInteraction(
      keyboard = list(enabled = TRUE),
      multiselect = TRUE,
      navigationButtons = FALSE
    ),
    configure = configure_component,
    legend = . %>% visNetwork::visLegend(
      enabled = TRUE,
      main = 'Legend',
      position = 'right',
      addEdges = edges_legend(network_representation_type),
      addNodes = nodes_legend(network_representation_type),
      useGroups = FALSE,
      zoom = TRUE
    )
  )
  return(network_components)
}
nodes_legend <- function(network_representation_type, config = setup_graph_config()) {
  return(
    data.frame(
      label = c("Celltype"),
      shape = c("dot"),
      color = c(config$NODE_COLORING$BACKGROUND))
  )
}
edges_legend <- function(network_representation_type, config = setup_graph_config()) {
  if (network_representation_type == 'ORA') {
    return(
      data.frame(
        color = c(
          config$EDGE_COLORING$COLOR_DOWN, config$EDGE_COLORING$COLOR_UP,
          config$EDGE_COLORING$COLOR_ROBUST, config$EDGE_COLORING$COLOR_DIFF
          ),
        label = c('DOWN', 'UP', 'STABLE', 'UP&DOWN'),
        arrows = c('to', 'to', 'to', 'to')
      )
    )
  
    } else if (network_representation_type %in% c("COUNTS_YOUNG", "COUNTS_OLD")) {
    return(
      data.frame(
        color = 'darkblue',
        label = '#interactions',
        arrows = 'to',
        width = 2
      )
    )
    } else if (network_representation_type == 'COUNTS_DIFF') {
      return(data.frame(
        color = 'darkblue',
        label = 'delta(interactions)',
        arrows = 'to',
        width = 2,
        length = 10
        )
      )
    } else {
      stop()
  }
}
build_network <- function(network_skeleton,
                          nodes_global = NULL,
                          edges_global = NULL,
                          layout = NULL,
                          legend = NULL,
                          options = NULL,
                          interactive = NULL,
                          configure = NULL) {
  vis_funcs <- list(
    nodes_global = nodes_global,
    edges_global = edges_global,
    layout = layout,
    legend = legend,
    options = options,
    interactive = interactive,
    configure = configure
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
      (vis_funcs$configure)
  )
}
interactive_from_igraph <- function(object,
                                    G,
                                    network_representation_type,
                                    network_layout_type,
                                    exclude_nonsign) {
  if (network_representation_type == 'ORA') {
    ne <- extract_node_edges(G, exclude_nonsign)
  } else if (network_representation_type %in% c('COUNTS_YOUNG', 'COUNTS_OLD', 'COUNTS_DIFF') ) {
    ne <- extract_node_edges(G, exclude_nonsign = FALSE)
  } else {
    stop('Invalid network_representation_type argument.')
  }
  nodes <- ne$nodes
  edges <- ne$edges
  
  network_components <- get_network_components(object, nodes, edges, G$layout, network_representation_type, network_layout_type)
  
  interactive_network <- do.call(build_network, network_components)
  return(interactive_network)
}
build_interactive_network <- function(object,
                                      network_representation_type = 'ORA',
                                      network_layout_type,
                                      class_signature,
                                      subobject_name
) {
  check_network_layout_type_arg(network_layout_type)
  G <- build_igraph(object, network_layout_type)

  interactive_net <- interactive_from_igraph(object, G, network_representation_type, network_layout_type, exclude_nonsign = TRUE)
  return(interactive_net)
}

setup_graph_config <- function() {
  # TODO: Review
  GRAPH_CONFIG <- list(
    EDGE_COLORING = list(
      COLOR_UP = "#277DA1", # blue
      COLOR_DOWN = "#F94144", # red
      COLOR_DIFF = "#F9C74F",
      COLOR_ROBUST = "#90BE6D",
      CUTOFF_UP = 1,  # TODO
      CUTOFF_DOWN = 1,
      CUTOFF_DIFF = 1,
      COLOR_NONE = rgb(0.2, 0.2, 0.2, alpha = 0.1)
    ),
    NODE_COLORING = list(
      BACKGROUND = '#f8961e',
      BORDER = '#577590',
      HIGHLIGHT = list(
        BACKGROUND = '#43AA8B',
        BORDER = '#577590'
      ),
      HOVER = list(
        BACKGROUND = '#43AA8B',
        BORDER = '#577590'
      )
    ),
    
    EDGE_STYLE = list(
      ARROW_SIZE = 0.5,
      WIDTH = 2.5
    ),

    LAYOUT = list(
      HGAP = 20,
      VGAP = 20
    ),

    VERTEX_STYLE = list(
      SIZE = 10,
      LABEL_DIST = 1.5,
      LABEL_CEX = 1.2,
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
    )
  )
  return(GRAPH_CONFIG)
}
