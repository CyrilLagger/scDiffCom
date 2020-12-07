# make_empty_scDiffCom <- function() {
#   return(new("scDiffCom"))
# }
# make_empty_igraph <- function() {
#   return(igraph::make_empty_graph())
# }
# make_empty_visNet <- function() {
#   return(visNetwork::visIgraph(make_empty_igraph()))
# }
# 
# setOldClass("igraph", make_empty_igraph())
# setOldClass("visNetwork", NULL)
# 
# setClass('InteractiveNetsBuilder',
#          slots = c(
#            INTERACTIONS = 'scDiffCom',
#            IGRAPH = "igraph",
#            VISNET = 'visNetwork',
#            CONFIG = 'list'
#          ),
#          prototype = list(
#            INTERACTIONS = make_empty_scDiffCom(),
#            IGRAPH = make_empty_igraph(),
#            VISNET = NULL,
#            CONFIG = list()
#          )
# )
# 
# setValidity('InteractiveNetsBuilder', function(object) {
#   # TODO
#   return(TRUE)
# })
# 
# 
# setGeneric('step01_scDiffCom_to_igraph',
#            function(x) standardGeneric('step01_scDiffCom_to_igraph'),
#            signature = 'x'
# )
# setMethod('step01_scDiffCom_to_igraph', "InteractiveNetsBuilder", 
#           function(x) {
#             build_igraph()
# })



construct_graph <- function(object) {

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
  COL_P_VALUE_UP <- "P_VALUE_UP"
  COL_P_VALUE_DOWN <- "P_VALUE_DOWN"
  
  COL_BH_P_VALUE_DIFF <- "BH_P_VALUE_DIFF"
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
    'COUNTS_TOTAL' := COUNTS_VALUE_REGULATED_DIFF + COUNTS_VALUE_REGULATED_FLAT
  ]
  
  subset_celltypes_dt_columns <- function(dt, cols) {
    if (length(cols) != 16) {
      stop('Not implemented')
    }
    dt <- dt[, list(
      get(cols[1]), get(cols[2]), get(cols[3]),
      get(cols[4]), get(cols[5]), get(cols[6]),
      get(cols[7]), get(cols[8]), get(cols[9]),
      get(cols[10]), get(cols[11]), get(cols[12]),
      get(cols[13]), get(cols[14]), get(cols[15]),
      get(cols[16])
    )]
    names(dt) = cols
    return(dt)
  }
  
  cols_to_select <- c(
    "Ligand_cell", "Receptor_cell", "OR_DIFF", "OR_UP", "OR_DOWN",
    COL_P_VALUE_DIFF, COL_P_VALUE_UP, COL_P_VALUE_DOWN,
    COL_BH_P_VALUE_DIFF, COL_BH_P_VALUE_UP, COL_BH_P_VALUE_DOWN,
    'COUNTS_UP', 'COUNTS_DOWN', 'COUNTS_DIFF', 'COUNTS_FLAT', 'COUNTS_TOTAL'
  )
  dt_ctypes = subset_celltypes_dt_columns(dt_ctypes, cols_to_select)
  return(dt_ctypes)
}
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
  has_inf = sum(dt_edge == Inf) > 0
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
        is_P_VALUE_UP <- BH_P_VALUE_UP < 0.05
        is_P_VALUE_DOWN <- BH_P_VALUE_DOWN < 0.05
      } else {
        is_P_VALUE_UP <- P_VALUE_UP < 0.05
        is_P_VALUE_DOWN <- P_VALUE_DOWN < 0.05
      }
      
      is_OR_UP_CHANGE <- OR_UP >= CUTOFF_CHANGE_UP
      is_OR_DOWN_CHANGE <- OR_DOWN >= CUTOFF_CHANGE_DOWN
      is_OR_UP_ROBUST <- OR_UP <= CUTOFF_ROBUST_UP
      is_OR_DOWN_ROBUST <- OR_DOWN <= CUTOFF_ROBUST_DOWN
      
      if (is_OR_UP_CHANGE & is_OR_DOWN_CHANGE & is_P_VALUE_UP & is_P_VALUE_DOWN) {
        return(LABEL_DIFF)
      } else if (is_OR_UP_CHANGE & is_P_VALUE_UP) {
        return(LABEL_UP)
      } else if (is_OR_DOWN_CHANGE & is_P_VALUE_DOWN) {
        return(LABEL_DOWN)
      } else if (is_OR_UP_ROBUST & is_OR_DOWN_ROBUST & is_P_VALUE_UP & is_P_VALUE_DOWN) {
        return(LABEL_ROBUST)
      } else {
        return(LABEL_NONE)
      }
    }
  )
}
map_edgetype_to_color <- function(edge_label, config = setup_graph_config()) {
  
  # TODO: Inline variables
  COLOR_DIFF <- config$EDGE_COLORING$COLOR_BOTH
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
infer_vertex_types <- function(vertex_names) {
  vertex_types <- purrr::map_lgl(
    strsplit(vertex_names, "[()]"),
    ~ .x[2] == "L"
  )
  return(vertex_types)
}
setup_graph <- function(G, config = setup_graph_config(), use_adjpval, disperse) {
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
add_vertex_size <- function(G) {
  MAXSIZE <- 20
  MINSIZE <- 5
  igraph::V(G)$vertex.size <- rescale(V(G)$counts, MINSIZE, MAXSIZE)
  return(G)
}
num_interactions_cellpairs <- function(object) {
  return(
    object@cci_detected[,
      .(count = .N),
      by = .(EMITTER_CELLTYPE, RECEIVER_CELLTYPE)
    ][, .(
      "Ligand_celltype" = EMITTER_CELLTYPE,  # TODO: Change these colnames
      "Receptor_celltype" = RECEIVER_CELLTYPE,
      "num_interacts" = count
    )]
  )
}
add_edge_width <- function(G) {
  WIDTH_MIN <- 0
  WIDTH_MAX <- 10
  num_interacts <- igraph::as_data_frame(G, 'edges')$COUNTS_TOTAL
  igraph::E(G)$width <- rescale(sqrt(num_interacts), WIDTH_MIN, WIDTH_MAX)
  return(G)
}
setup_vertices <- function(G, config) {
  igraph::V(G)$vertex_types <- infer_vertex_types(igraph::V(G)$name)
  G <- add_vertex_size(G)
  return(G)
}
setup_edges <- function(G, config, use_adjpval) {
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
rescale <- function(v, min_, max_) {
  (v - min(v)) / (max(v) - min(v)) * (max_ - min_) + min_
}

#### Functions for improved layout of biparite type ####
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
sort_layout_vertices <- function(layout, G, config, disperse) {
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
check_network_type_arg <- function(network_type) {
  if (!(network_type %in% c("bipartite", "celltypes"))) {
    message <- sprintf('network_type "%s" not in "bipartite", "celltypes"', network_type)
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
    vertex.counts = mean(vertex.counts)  # TODO: Not sure still works
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
build_igraph <- function(object,
                         network_type) {
  G <- construct_graph(object)
  G <- setup_graph(G, use_adjpval = TRUE, disperse = TRUE)
  if (network_type == "celltypes") {
    G <- map_bipartite_to_community(G)
  }
  return(G)
}
extract_node_edges <- function(G, exclude_nonsign = TRUE) {
  nodes <- data.table::data.table(igraph::as_data_frame(G, what = "vertices"))
  edges <- data.table::data.table(igraph::as_data_frame(G, what = "edges"))

  if (exclude_nonsign) {
    # TODO: Move to construct_graph
    edges <- edges[edge_type != "none"]
  }

  return(
    list(nodes = nodes, edges = edges)
  )
}
build_network_skeleton <- function(
                                 nodes,
                                 edges,
                                 network_type) {
  nodes_ <- data.table::copy(nodes)
  edges_ <- data.table::copy(edges)
  
  # TODO: Extract function
  node_annotation_html <- function() {
      temp_dt = data.table::data.table(
          Type=c("Tissue counts", "Cell counts", "Cell OR"),
          Total=1:3,
          UP=1:3,
          DOWN=1:3,
          FLAT=1:3,
          DIFF=1:3
      )
      table_html = kableExtra::kbl(temp_dt)
      html = paste0("<p>CELLNAME</p>", as.character(table_html), "<p>Num cells: ?</p>")
    return(html)
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
    "color" := "green"
  ][
    , # Add: color.background, color.border, color.highlight.background/border, color.hover.background/hover
    "shadow" := TRUE
  ][
    ,
    # 'shapeProperties' := NULL][,
    "title" := node_annotation_html()
  ]

  if (network_type == "bipartite") {
    nodes_[
      ,
      "group" := vertex_types
    ][
      ,
      "level" := ifelse(vertex_types == TRUE, 1, 2) # For hierarchical layout from vis
    ]
  } else if (network_type == "celltypes") {
  } else {
    stop('argument network_type not in "bipartite", "celltypes"')
  }
  
  edge_annotation_html <- function() {
    paste0("<p>E cell > R cell</p>", as.character(kableExtra::kbl(data.table::data.table(Type=c('Counts', 'OR'), Total=1:2, UP=1:2, DOWN=1:2, FLAT=1:2, DIFF=1:2))) )
  }
  edges_[
    ,
    color.color := edge.color
  ][
    ,
    color.highlight := edge.color
  ][
    ,
    color.hover := edge.color
  ][
    ,
    smooth := TRUE
  ][
    ,
    title := edge_annotation_html()  # "Num_UP=10, Num_DOWN=10, Num_FLAT=30"
  ]

  return(
    visNetwork::visNetwork(
      nodes_, edges_,
      width = "100%", # height='100%',
      main = sprintf("%s", 'tissue name'),
      submain = sprintf("%s", network_type),
      footer = "Ageing of intercellular communication",
      background = "white"
    )
  )
}
get_network_components <- function(
                                   nodes,
                                   edges,
                                   layout,
                                   network_type,
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
      nodes, edges,
      network_type = network_type
    ),
    nodes_global = . %>% visNetwork::visNodes(
      shape = "dot",
      physics = FALSE
    ),
    edges_global = . %>% visNetwork::visEdges(
      shadow = TRUE,
      arrows = "middle",
      smooth = list(enabled = TRUE, roundness = 0.5)
    ),
    layout = . %>% visNetwork::visIgraphLayout(
      layout = "layout.norm",
      layoutMatrix = layout
    ),
    options = . %>% visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      selectedBy = list(variable = "id", multiple = TRUE),
      manipulation = TRUE
    ),
    interactive = . %>% visNetwork::visInteraction(
      keyboard = list(enabled = TRUE),
      multiselect = TRUE,
      navigationButtons = TRUE
    ),
    configure = configure_component,
    legend = . %>% visNetwork::visLegend(
      enabled = TRUE,
      useGroups = FALSE
    )
  )
  return(network_components)
}
build_network <- function(
                          network_skeleton,
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
interactive_from_igraph <- function(
                                    G,
                                    network_type,
                                    exclude_nonsign = TRUE) {
  ne <- extract_node_edges(G)
  nodes <- ne$nodes
  edges <- ne$edges

  network_components <- get_network_components(nodes, edges, G$layout, network_type)
  interactive_network <- do.call(build_network, network_components)
  return(interactive_network)
}
build_interactive_network <- function(
  object,
  network_type,
  class_signature,
  subobject_name
) {
  check_network_type_arg(network_type)
  G <- build_igraph(object, network_type)

  interactive_net <- interactive_from_igraph(G, network_type)
  return(interactive_net)
}

setup_graph_config <- function() {
  # TODO: Review
  GRAPH_CONFIG <- list(
    EDGE_COLORING = list(
      COLOR_UP = "#4285F4", # blue
      CUTOFF_UP = 1.1,
      COLOR_DOWN = "#DB4437", # red
      CUTOFF_DOWN = 1.1,
      COLOR_BOTH = "#F4B400",
      COLOR_ROBUST = "#0F9D58",
      COLOR_NONE = rgb(0.2, 0.2, 0.2, alpha = 0.1)
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
