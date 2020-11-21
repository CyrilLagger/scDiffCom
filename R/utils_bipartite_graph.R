######

construct_graph <- function(
                            ora_ct,
                            cci_table_filtered,
                            graph_name) {
  dt_edge <- get_celltypes_enrichment(ora_ct)
  dt_edge <- dt_edge[, list(
    "Ligand_cell" = paste0(Ligand_cell, " (L)"),
    "Receptor_cell" = paste0(Receptor_cell, " (R)"),
    "OR_DIFF" = OR_DIFF,
    "OR_UP" = OR_UP,
    "OR_DOWN" = OR_DOWN,
    "pval_DIFF" = pval_DIFF,
    "pval_UP" = pval_UP,
    "pval_DOWN" = pval_DOWN,
    "pval_adj" = pval_readj_on_celltypes,
    "pval_adj_UP" =  pval_readj_on_celltypes_UP,
    "pval_adj_DOWN" = pval_readj_on_celltypes_DOWN
  )]
  dt_edge[is.na(dt_edge)] <- 1
  G <- igraph::graph_from_data_frame(
    d = dt_edge,
    directed = TRUE,
    vertices = NULL
  )
  G$graph_name <- graph_name
  G$ora <- ora_ct
  G$dt <- cci_table_filtered
  G$interacts_counts <- count_interactions_cellpairs_tissue(cci_table_filtered)
  G$celltype_counts <- count_celltypes_tissue(cci_table_filtered)
  return(G)
  # Subset the tissue
  # if( tissue == 'All' ) {
  #  stop('construct_graph: Not implemented for All. Gets in conflict with some
  #       downstream merging from filtered datasets that don\'t represent
  #       explicitely a tissue named "All"')
  # }
  # dt_edge = dt_ctypes[Tissue == tissue, .SD, .SDcols = !c('Tissue')]
  # message('construct_graph: colnames are hard-coded in get_celltypes_enrichment')
  # Label the transmitter (L) and receiver (R) cells
  ## Handle NAs
  # message(paste0('process_edge_dt: #NAs in ORA:celltypes for ',
  #               tissue, ':', sum(is.na(dt_edge)),'; NA->1'))
  # Handle Inf
  # Shouldn't have Inf with pseudocounts
  # CLIP_VAL = 50
  # if (sum(dt_edge==Inf) > 0) {
  #  stop('construct_graph: Inf values in dt_edge')
  # message(paste0('process_edge_dt: Inf value present in edge data.table for',
  #                tissue, '; Clipping to ', CLIP_VAL))
  # dt_edge[dt_edge == Inf] = min(dt_edge[, list(OR, OR_UP, OR_DOWN)], CLIP_VAL)
  # }
  # message('construct_graph: igraph G holds as graph attributes dt_ora and dt_filtered.')
}

get_celltypes_enrichment <- function(
                                     dt) {
  # message('get_celltype_enrichment: upgrade warning: OR -> OR_DIFF')
  COL_LIGAND_RECEPTOR_CELLTYPES <- "LR_CELLTYPE"
  COL_pval <- "pval_DIFF"
  COL_pval_UP <- "pval_UP"
  COL_pval_DOWN <- "pval_DOWN"
  COL_pval_readj <- "pval_readj_on_celltypes"
  COL_pval_readj_UP <- "pval_readj_on_celltypes_UP"
  COL_pval_readj_DOWN <- "pval_readj_on_celltypes_DOWN"
  dt_ctypes <- dt[
    Category == COL_LIGAND_RECEPTOR_CELLTYPES
  ][
    ,
    c(COL_pval_readj, COL_pval_readj_UP, COL_pval_readj_DOWN) := list(
      p.adjust(get(COL_pval), method = "BH"),
      p.adjust(get(COL_pval_UP), method = "BH"),
      p.adjust(get(COL_pval_DOWN), method = "BH")
    )
  ]
  dt_ctypes <- dt_ctypes[, c("Ligand_cell", "Receptor_cell") := list(
    sub("_.*", "", Value), # substitute w/ regex
    sub(".*_", "", Value)
  )]
  cols_to_select <- c(
    "Ligand_cell", "Receptor_cell", "OR_DIFF", "OR_UP", "OR_DOWN",
    COL_pval, COL_pval_UP, COL_pval_DOWN,
    COL_pval_readj, COL_pval_readj_UP, COL_pval_readj_DOWN
  )
  dt_ctypes <- dt_ctypes[, cols_to_select, with = FALSE]
  return(dt_ctypes)
}

count_interactions_cellpairs_tissue <- function(
                                                dt_filtered) {
  return(
    dt_filtered[
      ,
      .(count = .N),
      by = .(L_CELLTYPE, R_CELLTYPE)
    ][, .(
      "Ligand_celltype" = L_CELLTYPE,
      "Receptor_celltype" = R_CELLTYPE,
      "num_interacts" = count
    )]
  )
}

count_celltypes_tissue <- function(
                                   dt) {
  dt_t <- copy(dt)
  counts_dt <- funion(
    dt_t[, .(
      name = L_CELLTYPE,
      counts = (L_NCELLS_OLD + L_NCELLS_YOUNG) / 2 # average counts in young-old
    )],
    dt_t[, .(
      name = R_CELLTYPE,
      counts = (R_NCELLS_OLD + R_NCELLS_YOUNG) / 2 # average counts in young-old
    )]
  )
  counts_dt <- counts_dt[order(counts_dt$name)]
  return(counts_dt)
}

############

define_graph_config <- function() {
  GRAPH_CONFIG <- list(
    EDGE_COLORING = list(
      COLOR_UP = "#DB4437", # red
      CUTOFF_UP = 1.1,
      COLOR_DOWN = "#4285F4", # blue
      CUTOFF_DOWN = 1.1,
      COLOR_BOTH = "#F4B400",
      COLOR_ROBUST = "#0F9D58",
      # COLOR_NONE = "#dde2e4"
      COLOR_NONE = grDevices::rgb(0.2, 0.2, 0.2, alpha = 0.1)
    ),
    EDGE_STYLE = list(
      ARROW_SIZE = 0.5,
      WIDTH = 2.5 # will be adjusted by |OR|
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
}

############

setup_graph <- function(
                        G,
                        config,
                        use_adjpval,
                        disperse) {
  G <- setup_vertices(G, config)
  G <- setup_edges(G, config, use_adjpval)
  G <- setup_layout(G, config, disperse)
  return(G)
}

setup_vertices <- function(
                           G,
                           config) {
  igraph::V(G)$vertex_types <- infer_vertex_types(igraph::V(G)$name)
  G <- add_vertex_size(G)
  return(G)
}

infer_vertex_types <- function(
                               vertex_names) {
  # print(vertex_names)
  vertex_types <- purrr::map_lgl(
    strsplit(vertex_names, "[()]"),
    ~ .x[2] == "L"
  )
  # print(vertex_types)
  return(vertex_types)
}

add_vertex_size <- function(
                            G) {
  name <- counts <- NULL
  igraph::V(G)$vertex.size <- 30
  MAXSIZE <- 20
  MINSIZE <- 5
  vertex_names <- igraph::V(G)$name
  # counts_dt <- count_celltypes_tissue(
  #   G$dt,
  #   G$tissue
  # )
  celltypes_to_vertexnames <- function(counts_dt) {
    vertex_dt <- funion(
      counts_dt[, list(
        name = paste0(name, " (L)"),
        counts = counts
      )],
      counts_dt[, list(
        name = paste0(name, " (R)"),
        counts = counts
      )]
    )
    return(vertex_dt)
  }
  # vertex_counts_dt = celltypes_to_vertexnames(counts_dt)
  vertex_counts_dt <- celltypes_to_vertexnames(G$celltype_counts)
  vertex_counts_dt <- vertex_counts_dt[vertex_counts_dt$name %in% vertex_names]
  # Order data.table by vertex_names order and extract counts in correct order.
  ord <- match(vertex_names, vertex_counts_dt$name)
  vertex_counts <- vertex_counts_dt[ord, counts]
  igraph::V(G)$vertex.counts <- vertex_counts
  vertex_sizes <- scales::rescale(
    igraph::V(G)$vertex.counts,
    to = c(MINSIZE, MAXSIZE)
  )
  igraph::V(G)$vertex.size <- vertex_sizes
  return(G)
}

setup_edges <- function(
                        G,
                        config,
                        use_adjpval) {
  igraph::E(G)$arrow.size <- config$EDGE_STYLE$ARROW_SIZE
  igraph::E(G)$edge.color <- purrr::map_chr(igraph::E(G), ~ color_edge(.x, config = config, use_adjpval = use_adjpval))
  G <- add_edge_width(G)
  return(G)
}

color_edge <- function(
                       edge,
                       config = GRAPH_CONFIG,
                       use_adjpval = TRUE) {
  COLOR_UP_AND_DOWN <- config$EDGE_COLORING$COLOR_BOTH
  COLOR_UP <- config$EDGE_COLORING$COLOR_UP
  COLOR_DOWN <- config$EDGE_COLORING$COLOR_DOWN
  COLOR_ROBUST <- config$EDGE_COLORING$COLOR_ROBUST
  COLOR_NONE <- config$EDGE_COLORING$COLOR_NONE

  CUTOFF_CHANGE_UP <- config$EDGE_COLORING$CUTOFF_UP
  CUTOFF_CHANGE_DOWN <- config$EDGE_COLORING$CUTOFF_DOWN
  CUTOFF_ROBUST_UP <- 0.99
  CUTOFF_ROBUST_DOWN <- 0.99
  # message('color_edge: some hard-coded cutoff.')
  # message(paste0(
  #  'color_edge: refactor to set edge state in the ',
  #  'edge data table to be reused within further analysis.')
  # )
  if (use_adjpval) {
    is_pval_UP <- edge$pval_adj_UP < 0.05
    is_pval_DOWN <- edge$pval_adj_DOWN < 0.05
  } else {
    is_pval_UP <- edge$pval_UP < 0.05
    is_pval_DOWN <- edge$pval_DOWN < 0.05
  }
  is_OR_UP_CHANGE <- edge$OR_UP >= CUTOFF_CHANGE_UP
  is_OR_DOWN_CHANGE <- edge$OR_DOWN >= CUTOFF_CHANGE_DOWN

  is_OR_UP_ROBUST <- edge$OR_UP <= CUTOFF_ROBUST_UP
  is_OR_DOWN_ROBUST <- edge$OR_DOWN <= CUTOFF_ROBUST_DOWN

  if (is_OR_UP_CHANGE & is_OR_DOWN_CHANGE & is_pval_UP & is_pval_DOWN) {
    return(COLOR_UP_AND_DOWN)
  } else if (is_OR_UP_CHANGE & is_pval_UP) {
    return(COLOR_UP)
  } else if (is_OR_DOWN_CHANGE & is_pval_DOWN) {
    return(COLOR_DOWN)
  } else if (is_OR_UP_ROBUST & is_OR_DOWN_ROBUST & is_pval_UP & is_pval_DOWN) {
    return(COLOR_ROBUST)
  } else {
    return(COLOR_NONE)
  }

  # if ((edge$OR_UP >= CUTOFF_UP) & (edge$OR_DOWN >= CUTOFF_DOWN)) {
  #   return(COLOR_UP_AND_DOWN)
  # } else if (edge$OR_UP >= CUTOFF_UP) {
  #   return(COLOR_UP)
  # } else if (edge$OR_DOWN >= CUTOFF_DOWN) {
  #   return(COLOR_DOWN)
  # } else {
  #   return(COLOR_NONE)
  # }
}

add_edge_width <- function(
                           G) {
  Ligand_celltype <- Receptor_celltype <- NULL
  WIDTH_MIN <- 0
  WIDTH_MAX <- 10
  # or_max = mapply(max, igraph::E(G)$OR_UP, igraph::E(G)$OR_DOWN)
  edge_dt <- data.table(
    igraph::as_data_frame(G)[, c("from", "to")]
  )
  names(edge_dt) <- c("Ligand_celltype", "Receptor_celltype")
  edge_dt[, Ligand_celltype :=
    sub("(*) [(]L[)]", "", Ligand_celltype)][
    ,
    Receptor_celltype :=
      sub("(*) [(]R[)]", "", Receptor_celltype)
  ]
  edge_dt <- merge(
    edge_dt, G$interacts_counts,
    by = c("Ligand_celltype", "Receptor_celltype"),
    all.x = TRUE,
    sort = FALSE
  )
  edge_dt[is.na(edge_dt)] <- 0
  # message('add_edge_width: Some NAs occur, unclear reason. For now NA->0.')
  # Correct order is accomplished by the merge.
  num_interacts <- edge_dt[, num_interacts]
  igraph::E(G)$width <- scales::rescale(
    sqrt(num_interacts),
    to = c(WIDTH_MIN, WIDTH_MAX)
  )
  # message('add_edge_width: widths are scaled non-linearly.')
  return(G)
}

setup_layout <- function(
                         G,
                         config,
                         disperse) {
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


#############

plot_graph <- function(
                       G,
                       config = GRAPH_CONFIG,
                       path = NULL,
                       show_legend = FALSE) {
  MAIN <- paste0("Communication profile with age: ", G$graph_name)
  MARGIN <- 0.5
  SUBTITLE <- ""
  # SUBTITLE = paste0('Difference graph of overrepresented celltype ',
  #                   'communication that get\n altered with ageing')
  if (!is.null(path)) {
    pdf(path)
  }
  igraph::plot.igraph(
    G,
    layout = G$layout,
    vertex.color = config$VERTEX_STYLE$COLOR,
    vertex.size = igraph::V(G)$vertex.size,
    vertex.label.dist = config$VERTEX_STYLE$LABEL_DIST,
    vertex.label.cex = config$VERTEX_STYLE$LABEL_CEX,
    edge.color = igraph::E(G)$edge.color,
    edge.width = igraph::E(G)$edge.width,
    main = MAIN,
    margin = MARGIN,
    sub = SUBTITLE
  )
  if (show_legend) {
    LEGEND_TITLE <- "Edge color legend"
    LEGEND_COLORS <- c(
      config$EDGE_COLORING$COLOR_NONE,
      config$EDGE_COLORING$COLOR_UP,
      config$EDGE_COLORING$COLOR_DOWN,
      config$EDGE_COLORING$COLOR_BOTH
    )

    legend(
      x = -1.5, y = -1.1,

      legend = config$LEGEND$LEGEND_LABELS,
      col = LEGEND_COLORS,
      title = LEGEND_TITLE,

      bty = "o",
      pch = config$LEGEND$PCH,
      cex = config$LEGEND$CEX,
      pt.cex = config$LEGEND$PT.CEX,
      bg = config$LEGEND$BG,
      ncol = config$LEGEND$NCOL
    )
  }
  if (!is.null(path)) {
    dev.off()
  }
}


create_analysis_dirs <- function(
                                 dir,
                                 analysis_name,
                                 subdirs) {
  analysis_dir <- file.path(dir, analysis_name)
  dir.create(analysis_dir)

  map(
    subdirs,
    ~ dir.create(file.path(analysis_dir, .x))
  )
}

weight_edge <- function(
                        edge,
                        min_width,
                        max_width,
                        min_value,
                        max_value,
                        config = GRAPH_CONFIG) {

  # TODO: weight_edge and color edge can be combined in process_edge
  # TODO: scale edge width by OR magnitude?
  # (c + min_width)*(max_width/min_width)

  CUTOFF_UP <- config$EDGE_COLORING$CUTOFF_UP
  CUTOFF_DOWN <- config$EDGE_COLORING$CUTOFF_DOWN

  # CUTOFF_CHANGE_UP = config$EDGE_COLORING$CUTOFF_UP
  # CUTOFF_CHANGE_DOWN = config$EDGE_COLORING$CUTOFF_DOWN
  # CUTOFF_ROBUST_UP = 0.99
  # CUTOFF_ROBUST_DOWN = 0.99
  # message('weight_edge: some hard-coded cutoff.')
  #
  # is_pval_UP = edge$pval_adj_UP < 0.05
  # is_pval_DOWN = edge$pval_adj_DOWN < 0.05
  #
  # is_OR_UP_CHANGE = edge$OR_UP >= CUTOFF_CHANGE_UP
  # is_OR_DOWN_CHANGE = edge$OR_DOWN >= CUTOFF_CHANGE_DOWN
  #
  # is_OR_UP_ROBUST = edge$OR_UP <= CUTOFF_ROBUST_UP
  # is_OR_DOWN_ROBUST = edge$OR_DOWN <= CUTOFF_ROBUST_DOWN


  if ((edge$OR_UP >= CUTOFF_UP) & (edge$OR_DOWN >= CUTOFF_DOWN)) {
    return(edge$OR)
  } else if (edge$OR_UP >= CUTOFF_UP) {
    return(edge$OR_UP)
  } else if (edge$OR_DOWN >= CUTOFF_DOWN) {
    return(edge$OR_DOWN)
  } else {
    return(1)
  }
}


get_edges_from_vertex <- function(
                                  v_name,
                                  G) {
  # return(igraph::E(G)[from(v_name)])
  return(igraph::incident(G, v_name, mode = "out"))
}

get_edges_to_vertex <- function(
                                v_name,
                                G) {
  return(igraph::incident(G, v_name, mode = "in"))
}

count_up_and_down_edges <- function(
                                    edges,
                                    config) {
  color_up <- config$EDGE_COLORING$COLOR_UP
  color_down <- config$EDGE_COLORING$COLOR_DOWN

  counts <- table(edges$edge.color)
  num_up <- as.integer(counts[color_up])
  num_down <- as.integer(counts[color_down])

  num_up <- ifelse(is.na(num_up), 0, num_up)
  num_down <- ifelse(is.na(num_down), 0, num_down)

  return(list(N_UP = num_up, N_DOWN = num_down))
}

vertex_sort_key_atomic <- function(
                                   v_name,
                                   from,
                                   G,
                                   config) {
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

vertex_sort_keys <- function(
                             v_names,
                             from,
                             G,
                             config) {
  sort_keys <- purrr::map_dbl(
    v_names,
    ~ vertex_sort_key_atomic(.x, from, G, config)
  )
  return(sort_keys)
}

sort_layout_vertices <- function(
                                 layout,
                                 G,
                                 config,
                                 disperse) {
  vertex_names <- igraph::V(G)$name
  num_v <- length(vertex_names)
  midpoint <- num_v / 2
  ligand_vertices <- vertex_names[1:midpoint]
  receptor_vertices <- vertex_names[(midpoint + 1):num_v]

  ligand_keys <- vertex_sort_keys(ligand_vertices, from = TRUE, G, config)
  receptor_keys <- vertex_sort_keys(receptor_vertices, from = FALSE, G, config)

  # TODO: refactor into clear code
  # The reorders layout rows to achieve sorting
  layout_copy <- copy(layout)
  new_layout <- copy(layout)
  new_layout[order(ligand_keys), ] <- layout_copy[1:midpoint, ]
  new_layout[midpoint + order(receptor_keys), ] <- layout_copy[(midpoint + 1):num_v, ]

  if (disperse) {
    ## Make as function
    # Get hgap
    vgap <- config$LAYOUT$HGAP
    scale_factor <- 3

    new_layout[1:midpoint, 2][ligand_keys == 0] <- new_layout[1:9, 2][ligand_keys == 0] + scale_factor * vgap
    new_layout[1:midpoint, 2][ligand_keys > 0] <- new_layout[1:9, 2][ligand_keys > 0] + 2 * scale_factor * vgap

    new_layout[(midpoint + 1):num_v, 2][receptor_keys == 0] <- new_layout[(midpoint + 1):num_v, 2][receptor_keys == 0] + scale_factor * vgap
    new_layout[(midpoint + 1):num_v, 2][receptor_keys > 0] <- new_layout[(midpoint + 1):num_v, 2][receptor_keys > 0] + 2 * scale_factor * vgap
  }
  # new_layout = disperse_layout_vertices(new_layout, G, config)
  return(new_layout)
}

write_as_edge_table <- function(
                                G,
                                path) {
  if (!is.null(path)) {
    df_edge <- igraph::as_data_frame(G, what = "edges")
    write.table(
      df_edge,
      file = path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  } else {
    stop("write_as_edge_table: path argument not provided.")
  }
}
